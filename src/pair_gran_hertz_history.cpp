/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "pair_gran_hertz_history.h"
#include "atom.h"
#include "force.h"
#include "neigh_list.h"
#include "error.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHertzHistory::PairGranHertzHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp)
{
  no_virial_compute = 1;
  history = 1;
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */

inline void PairGranHertzHistory::deriveContactModelParams(int &ip, int &jp,double &meff, double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu) 
{
    #define LMP_GRAN_DEFS_DEFINE
    #include "pair_gran_defs.h"
    #undef LMP_GRAN_DEFS_DEFINE

    double reff=ri*rj/(ri+rj);

    double sqrtval = sqrt(reff*deltan);

    double Sn=2.*mpg->Yeff[itype][jtype]*sqrtval;
    double St=8.*mpg->Geff[itype][jtype]*sqrtval;

    kn=4./3.*mpg->Yeff[itype][jtype]*sqrtval;
    kt=St;
    gamman=-2.*sqrtFiveOverSix*mpg->betaeff[itype][jtype]*sqrt(Sn*meff);
    gammat=-2.*sqrtFiveOverSix*mpg->betaeff[itype][jtype]*sqrt(St*meff);
    xmu=mpg->coeffFrict[itype][jtype];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    /* NP
    char *testmsg=new char[200];
    sprintf(testmsg,"Yeff=%f,reff=%f,deltan=%f, kn=%f, kt=%f, gamman=%f, gammat=%f, xmu=%f\n",Yeff, reff,deltan, kn,kt,gamman,gammat,xmu);
    error->warning(testmsg);
    delete []testmsg;*/
    #define LMP_GRAN_DEFS_UNDEFINE
    #include "pair_gran_defs.h"
    #undef LMP_GRAN_DEFS_UNDEFINE

    return;
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHertzHistory::init_substyle()
{
  mpg->getMaterialParams(0,cohesionflag);
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairGranHertzHistory::settings(int narg, char **arg) 
{
  if (narg != 2) error->all("Illegal pair_style command");

  dampflag = force->inumeric(arg[0]);
  cohesionflag = force->inumeric(arg[1]);

  if (dampflag < 0 || dampflag > 1 || cohesionflag < 0 || cohesionflag > 1)
    error->all("Illegal pair_style command");
}

void PairGranHertzHistory::emit_results(int step, char *fname) {
  int inum = list->inum;
  int nall = atom->nlocal + atom->nghost;

  FILE *ofile = fopen(fname, "a");
  fprintf(ofile, "Step %d\n", step);
  fprintf(ofile, "Force comparison\n");
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "f[%d][%d], %.16f\n", i,j,atom->f[i][j]);
    }
  }
  fprintf(ofile, "Torque comparison\n");
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "torque[%d][%d], %.16f\n", i,j,atom->torque[i][j]);
    }
  }
  fprintf(ofile, "Shear comparison\n");
  for (int ii=0; ii<inum; ii++) {
    int i = list->ilist[ii];
    int jnum = list->numneigh[i];

    for (int jj = 0; jj<jnum; jj++) {
      int j = list->firstneigh[i][jj];

      double *shear = &(listgranhistory->firstdouble[i][3*jj]);
      for (int k=0; k<3; k++) {
        fprintf(ofile, "shear[%d][%d][%d], %.16f\n", i,j,k,shear[k]);
      }
    }
  }
  fflush(ofile);
  fclose(ofile);
}

void PairGranHertzHistory::emit_particle_details(int i, bool do_header=true) {
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

void PairGranHertzHistory::compute(int eflag, int vflag) {
  static int step = -1; step++;
  if (step > 26) {
    printf("[C] Test run for step %d!\n", step);
  }

  PairGranHookeHistory::compute(eflag, vflag);

  if (step > 26) {
    emit_results(step, "cpu.out");
  }
  if (step == 100) exit(0);
}
