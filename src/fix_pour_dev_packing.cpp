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

#include "limits.h"
#include "math.h"
#include "stdlib.h"
#include "domain.h"
#include "error.h"
#include "fix_pour_dev_packing.h"
#include "random_park.h"

#define INFINITE static_cast<int>(LONG_MAX/2)

using namespace LAMMPS_NS;

FixPourDevPacking::FixPourDevPacking(LAMMPS *lmp, int narg, char **arg) :
  FixPourDev(lmp, narg, arg)
{
  if(massflowrate>0.) error->all("Illegal fix pour/dev/packing command, keyword 'massflowrate' can not be used for packing");
  if(rate!=0.) error->all("Illegal fix pour/dev/packing command, keyword 'rate' can not be used for packing");

  //only do one insertion; add infinite
  nfreq=INFINITE;
  ninsert=nper;

  //print stats
  if (me == 0) {
    if (screen)
      if(particles_per_insertion()==1) fprintf(screen,
	      "Particle packing insertion: trying %f particles once\n",
	      nper);
	  else fprintf(screen,
	      "Rigid body packing insertion: trying %f particles once\nEach rigid body consists of %d particles\n",
	      nper,particles_per_insertion());
    if (logfile)
      if(particles_per_insertion()==1) fprintf(logfile,
	      "Particle packing insertion: trying %f particles once\n",
	      nper);
	  else fprintf(logfile,
	      "Rigid body packing insertion: trying %f particles once\nEach rigid body consists of %d particles\n",
	      nper,particles_per_insertion());
  }

}

/* ---------------------------------------------------------------------- */

FixPourDevPacking::~FixPourDevPacking()
{
    //nothing to do here
}

/* ---------------------------------------------------------------------- */

void FixPourDevPacking::init()
{
  if (domain->triclinic) error->all("Cannot use fix pour with triclinic box");
  init_substyle();
}

/* ----------------------------------------------------------------------
   random insertion heights of particle
------------------------------------------------------------------------- */

inline void FixPourDevPacking::random_insert_height(double &h,double &tmp,double vzrel)
{
   h=(lo_current+shift_randompos())+random->uniform()*(hi_current-lo_current-2.*shift_randompos());
}

/* ----------------------------------------------------------------------
   calculate inlet velocity of particle
------------------------------------------------------------------------- */

inline void FixPourDevPacking::calc_insert_velocities(int i,double **xnear,double &vxtmp,double &vytmp,double &vztmp)
{
   if (domain->dimension == 3) {
      vxtmp = rand_pour(vxlo,vxhi,vel_rand_style); 
      vytmp = rand_pour(vylo,vyhi,vel_rand_style); 
      vztmp = - sqrt(vz*vz); 
    } else {
      vxtmp = rand_pour(vxlo,vxhi,vel_rand_style); 
      vytmp = vy - sqrt(2.0*grav*(hi_current-xnear[i][1]));  
      vztmp = 0.0;
    }
}
