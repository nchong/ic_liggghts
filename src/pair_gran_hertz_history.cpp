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

#include "hertz.pb.h"
void PairGranHertzHistory::emit(const char *fname) {
  int inum = list->inum;
  int nall = atom->nlocal + atom->nghost;

  FILE *ofile = fopen(fname, "a");
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "x[%d][%d], %.16f\n", i,j,atom->x[i][j]);
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "v[%d][%d], %.16f\n", i,j,atom->v[i][j]);
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "omega[%d][%d], %.16f\n", i,j,atom->omega[i][j]);
    }
  }
  for (int i=0; i<nall; i++) {
    fprintf(ofile, "radius[%d], %.16f\n", i,atom->radius[i]);
  }
  for (int i=0; i<nall; i++) {
    fprintf(ofile, "rmass[%d], %.16f\n", i,atom->rmass[i]);
  }
  for (int i=0; i<nall; i++) {
    fprintf(ofile, "type[%d], %d\n", i,atom->type[i]);
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "f[%d][%d], %.16f\n", i,j,atom->f[i][j]);
    }
  }
  for (int i=0; i<nall; i++) {
    for (int j=0; j<3; j++) {
      fprintf(ofile, "torque[%d][%d], %.16f\n", i,j,atom->torque[i][j]);
    }
  }
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

#include "simple_timer.h"
static SimpleTimer *timers = new SimpleTimer[8];

void PairGranHertzHistory::compute(int eflag, int vflag) {
  static int step = -1; step++;
  static int cpu_steps = 0;

  if (step == 1000) {
    emit("step1000.in");
    PairGranHookeHistory::compute(eflag, vflag);
    emit("step1000.out");
    exit(0);
  }

  int inum = list->inum;
  if (inum != 0) { cpu_steps++; }
  timers[0].start();
  PairGranHookeHistory::compute(eflag, vflag);
  timers[0].stop();
  timers[0].add_to_total();
  if (cpu_steps % 1000 == 0) {
    printf("- [total time] %.3fms\n", timers[0].total_time());
    timers[0].reset();
  }
}
