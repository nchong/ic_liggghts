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

#include <string>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// External functions from cuda library for force decomposition

bool hertz_gpu_init(double *boxlo, double *boxhi, double cell_len, double skin);
void hertz_gpu_clear();
double hertz_gpu_cell();
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

