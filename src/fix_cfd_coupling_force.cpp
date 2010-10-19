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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling_force.h"
#include "fix_propertyPerAtom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::FixCfdCouplingForce(LAMMPS *lmp, int narg, char **arg) :  FixCfdCoupling(lmp, narg, arg)
{
    dragforce = NULL;
}

FixCfdCouplingForce::~FixCfdCouplingForce()
{
    if(dragforce) modify->delete_fix("dragforce");

}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingForce::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::init_submodel()
{
  special_settings();
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::special_settings()
{
  //register dragforce
  if(!dragforce)
  {
        char* fixarg[11];
        fixarg[0]="dragforce";
        fixarg[1]="all";
        fixarg[2]="property/peratom";
        fixarg[3]="dragforce";
        fixarg[4]="vector"; 
        fixarg[5]="no";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        dragforce = modify->add_fix_property_peratom(11,fixarg);
  }

  //values to be transfered to OF
  add_push_property("x","vector");
  add_push_property("v","vector");
  add_push_property("radius","scalar");

  //should make error check on push properties here
  
  //values to come from OF
  add_pull_property("dragforce","vector");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **dragf = dragforce->array_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      f[i][0] += dragf[i][0];
      f[i][1] += dragf[i][1];
      f[i][2] += dragf[i][2];
      
    }
  
}
