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
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling_force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::FixCfdCouplingForce(LAMMPS *lmp, int narg, char **arg) :  FixCfdCoupling(lmp, narg, arg)
{
      special_settings();
}

FixCfdCouplingForce::~FixCfdCouplingForce()
{

}

void FixCfdCouplingForce::special_settings()
{
  //values to be transfered to OF
  add_push_property("x","vector");
  add_push_property("v","vector");
  add_push_property("radius","scalar");

  //should make error check on push properties here
  
  //values to come from OF - will be stored in array_atom, beginning from 0
  dragforce_index = add_pull_property("dragforce","vector");

  size_peratom_cols = nvalues;

  grow_arrays(atom->nmax);

  zero_arrays();
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

void FixCfdCouplingForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      f[i][0] += array_atom[i][dragforce_index+0];
      f[i][1] += array_atom[i][dragforce_index+1];
      f[i][2] += array_atom[i][dragforce_index+2];
      
    }
  
}
