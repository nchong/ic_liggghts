/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "assert.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_gran_hertz_history_gpu.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "atom_vec.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"
#include "universe.h"
#include "mech_param_gran.h"

#include <string>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// External functions from cuda library for force decomposition

bool hertz_gpu_init(double *boxlo, double *boxhi, double cell_len, double skin);
void hertz_gpu_clear();
struct hashmap **create_shearmap(
  int inum, /*list->inum*/
  int *ilist, /*list->ilist*/
  int *numneigh, /*list->numneigh*/
  int **firstneigh, /*list->firstneigh*/
  int **firsttouch, /*listgranhistory->firstneigh*/
  double **firstshear /*listgranhistory->firstdouble*/);
void compare_shearmap(
  FILE *ofile,
  struct hashmap **shearmap,
  int inum, /*list->inum*/
  int *ilist, /*list->ilist*/
  int *numneigh, /*list->numneigh*/
  int **firstneigh, /*list->firstneigh*/
  int **firsttouch, /*listgranhistory->firstneigh*/
  double **firstshear /*listgranhistory->firstdouble*/);
double hertz_gpu_cell(
  const bool eflag, const bool vflag,
  const int inum, const int nall, const int ago,

  //inputs
  double **host_x, double **host_v, double **host_omega,
  double *host_radius, double *host_rmass,
  int *host_type,
  double dt,

  //stiffness parameters
  int num_atom_types,
  double **host_Yeff,
  double **host_Geff,
  double **host_betaeff,
  double **host_coeffFrict,
  double nktv2p,

  //inouts 
  struct hashmap **&host_shear, double **host_torque, double **host_force);
void hertz_gpu_name(const int gpu_id, const int max_nbors, char * name);
void hertz_gpu_time();
double hertz_gpu_bytes();

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHertzHistoryGPU::PairGranHertzHistoryGPU(LAMMPS *lmp) : PairGranHertzHistory(lmp)
{
}

PairGranHertzHistoryGPU::~PairGranHertzHistoryGPU()
{
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranHertzHistoryGPU::init_style()
{
  PairGranHertzHistory::init_style();

  // set the GPU ID
  int my_gpu=comm->me;
  my_gpu=0;
  int ngpus=1;

  // Repeat cutsq calculation because done after call to init_style
  double cut;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }

  // use the max cutoff length as the cell length
  double maxcut = -1.0;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if (cutsq[i][j] > maxcut) maxcut = cutsq[i][j];

  // for this problem, adding skin results in better perf
  // this may be a parameter in the future
  double cell_len = sqrt(maxcut) + neighbor->skin;

  if (!hertz_gpu_init(domain->boxlo, domain->boxhi, cell_len, neighbor->skin))
    error->one("At least one process could not allocate a CUDA-enabled gpu");
  
  if (comm->me == 0 && screen) {
    printf("\n-------------------------------------");
    printf("-------------------------------------\n");
    printf("- Using GPGPU acceleration for Hertz:\n");
    printf("-------------------------------------");
    printf("-------------------------------------\n");

    for (int i=0; i<ngpus; i++) {
      int gpui=my_gpu;
      gpui=i;
      char gpu_string[500];
      hertz_gpu_name(gpui,neighbor->oneatom,gpu_string);
      printf("GPU %d: %s\n",gpui,gpu_string);   
    }
    printf("-------------------------------------");
    printf("-------------------------------------\n\n");
  }
}

void PairGranHertzHistoryGPU::emit_particle_details(int i, bool do_header=true) {
  if (do_header) {
    printf("particle %d\n", i);
    printf("x\n%.16f\n%.16f\n%.16f\n",
      atom->x[i][0], atom->x[i][1], atom->x[i][2]);
    printf("v\n%.16f\n%.16f\n%.16f\n",
      atom->v[i][0], atom->v[i][1], atom->v[i][2]);
    printf("omega\n%.16f\n%.16f\n%.16f\n",
      atom->omega[i][0], atom->omega[i][1], atom->omega[i][2]);
    printf("radius\n%.16f\n", atom->radius[i]);
    printf("rmass\n%.16f\n", atom->rmass[i]);
  }

  printf("torque\n%.16f\n%.16f\n%.16f\n",
    atom->torque[i][0], atom->torque[i][1], atom->torque[i][2]);
  printf("force\n%.16f\n%.16f\n%.16f\n",
    atom->f[i][0], atom->f[i][1], atom->f[i][2]);
}

void PairGranHertzHistoryGPU::compute(int eflag, int vflag) {
  int inum = list->inum;
  static int step = -1; step++;

  if (inum == 0 || step < 27) {
    PairGranHertzHistory::compute(eflag, vflag);
  } else {
    printf("Test run for step %d!\n", step);

    inum = atom->nlocal;
    int nall = atom->nlocal + atom->nghost;
    int num_atom_types =  mpg->max_type+1;

    //inouts: shear, torque and force
    //create deep copies so we can compare results with CPU simulation
    struct hashmap **shearmap = create_shearmap(
        inum, list->ilist, list->numneigh, list->firstneigh,
        listgranhistory->firstneigh, listgranhistory->firstdouble);

    size_t table_size = nall * 3 * sizeof(double);
    double **gpu_torque = (double **)malloc(table_size);
    memcpy(gpu_torque, atom->torque, table_size);
    double **gpu_force = (double **)malloc(table_size);
    memcpy(gpu_force, atom->f, table_size);

    //paranoid
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        if (gpu_torque[i][j] != atom->torque[i][j]) {
          printf("torque[%d][%d] mismatch: %f / %f\n", 
              i, j, gpu_torque[i][j], atom->torque[i][j]);
          exit(1);
        }
        if (gpu_force[i][j] != atom->f[i][j]) {
          printf("force[%d][%d] mismatch: %f / %f\n", 
              i, j, gpu_force[i][j], atom->f[i][j]);
          exit(1);
        }
      }
    }

    //call gpu into dummy results (shear, torque and shear)
    hertz_gpu_cell(eflag, vflag, inum, nall, neighbor->ago,
      atom->x, atom->v, atom->omega,
      atom->radius, atom->rmass,
      atom->type,
      dt,

      num_atom_types,
      mpg->Yeff,
      mpg->Geff,
      mpg->betaeff,
      mpg->coeffFrict,
      force->nktv2p,

      shearmap,
      gpu_torque,
      gpu_force);

    //call cpu golden reference
    PairGranHertzHistory::compute(eflag, vflag);

    //printf("Post-kernel\n");
    FILE *ofile = fopen("compare.csv", "a");
    fprintf(ofile, "Step %d\n", step);
    fprintf(ofile, "Force comparison\n");
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        fprintf(ofile, "f[%d][%d], %.16f, %.16f\n", i,j,atom->f[i][j],gpu_force[i][j]);
      }
    }
    fprintf(ofile, "Torque comparison\n");
    for (int i=0; i<nall; i++) {
      for (int j=0; j<3; j++) {
        fprintf(ofile, "torque[%d][%d], %.16f, %.16f\n", i,j,atom->torque[i][j],gpu_torque[i][j]);
      }
    }
    fprintf(ofile, "Shear comparison\n");
    compare_shearmap(
        ofile, shearmap,
        inum, list->ilist, list->numneigh, list->firstneigh,
        listgranhistory->firstneigh, listgranhistory->firstdouble);

    fflush(ofile);
    fclose(ofile);

    exit(0);
    if (step > 30) {
      printf("The end, for now!\n");
      exit(0);
    }
  }
}
