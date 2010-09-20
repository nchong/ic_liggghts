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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "fix_propertyPerAtom.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001
#define myAtof lmp->force->numeric 

/* ---------------------------------------------------------------------- */

FixPropertyPerAtom::FixPropertyPerAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    //Check args
    if (narg <7) error->all("Illegal fix property/peratom command, not enough arguments");
    if (narg >17) error->warning("Vector length in fix property/peratom larger than 10. Are you sure you want that?");

    //Read args
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    if (strcmp(arg[4],"scalar") == 0) vectorStyle = 0;
    else if (strcmp(arg[4],"vector") == 0) vectorStyle = 1;
    else error->all("Unknown style for fix property/peratom. Valid styles are scalar or vector");

    if (strcmp(arg[5],"yes") == 0) restart_peratom = 1; //fix handles restart properties per atom
    else if (strcmp(arg[5],"no") == 0) restart_peratom = 0; //fix does not handle restart properties per atom
    else error->all("Unknown restart style for fix property/peratom. Valid styles are yes or no");

    if (strcmp(arg[6],"yes") == 0) commGhost = 1; //fix handles restart properties per atom
    else if (strcmp(arg[6],"no") == 0) commGhost = 0; //fix does not handle restart properties per atom
    else error->all("Unknown communicate_ghost style for fix property/peratom. Valid styles are yes or no");

    if (strcmp(arg[7],"yes") == 0) commGhostRev = 1; //fix handles restart properties per atom
    else if (strcmp(arg[7],"no") == 0) commGhostRev = 0; //fix does not handle restart properties per atom
    else error->all("Unknown communicate_reverse_ghost style for fix property/peratom. Valid styles are yes or no");

    nvalues = narg - 8;
    if ((nvalues==1) && (vectorStyle))
      error->all("Error in fix property/peratom: Number of default values provided not consistent with vector style. Provide more than 1 value or use style 'scalar'");
    defaultvalues = new double[nvalues];

    for (int j=0;j<nvalues;j++) defaultvalues[j] = myAtof(arg[8+j]);

    if (vectorStyle) size_peratom_cols = nvalues;
    else size_peratom_cols = 0;

    peratom_flag=1; 
    peratom_freq=1;
    extvector=0; 
    create_attribute = 1; 
    
    time_depend = 1; 

    if (commGhost) comm_forward = nvalues;
    if (commGhostRev) comm_reverse = nvalues;

    // perform initial allocation of atom-based array
    // register with Atom class
    vector_atom = NULL; array_atom = NULL;
    grow_arrays(atom->nmax); 
    atom->add_callback(0); 
    if (restart_peratom) atom->add_callback(1); 

    //zero all arrays since dump may access it on timestep 0
    //or a variable may access it before first run
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++){
      if (vectorStyle) for (int m = 0; m < nvalues; m++) array_atom[i][m] = 0.;
      else vector_atom[i]=0.;
    }

    //check if there is already a fix that tries to register a property with the same name
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyPerAtom*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->all("Error in fix property/peratom. There is already a fix that registers a variable of the same name");

}

/* ---------------------------------------------------------------------- */

FixPropertyPerAtom::~FixPropertyPerAtom()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  if (restart_peratom) atom->delete_callback(id,1);

  // delete locally stored arrays
  delete[] variablename;
  delete[] defaultvalues;

  if (vectorStyle) memory->destroy_2d_double_array(array_atom);
  else delete[] vector_atom;
}

/* ---------------------------------------------------------------------- */

int FixPropertyPerAtom::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ----------------------------------------------------------------------
   forward and backward comm to be used by other fixes as needed
------------------------------------------------------------------------- */

void FixPropertyPerAtom::do_forward_comm()
{
    if (commGhost) comm->forward_comm_fix(this);
    else error->all("forward_comm invoked, but not registered");
}

void FixPropertyPerAtom::do_reverse_comm()
{
   if (commGhostRev)  comm->reverse_comm_fix(this);
   else error->all("reverse_comm invoked, but not registered");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyPerAtom::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyPerAtom::grow_arrays(int nmax)
{
  if (vectorStyle) array_atom = memory->grow_2d_double_array(array_atom,nmax,nvalues,"FixPropertyPerAtom:array_atom");
  else vector_atom = (double*)(lmp->memory->srealloc(vector_atom, nmax*sizeof(double), "FixPropertyPerAtom:vector_atom"));
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyPerAtom::copy_arrays(int i, int j)
{
    if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[j][k]=array_atom[i][k];
    else vector_atom[j]=vector_atom[i];
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyPerAtom::set_arrays(int i)
{
    if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[i][k]=defaultvalues[k];
    else vector_atom[i]=defaultvalues[0];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyPerAtom::pack_exchange(int i, double *buf)
{
    if (vectorStyle) for(int k=0;k<nvalues;k++) buf[k] = array_atom[i][k];
    else buf[0] = vector_atom[i];
    return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixPropertyPerAtom::unpack_exchange(int nlocal, double *buf)
{
    if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[nlocal][k]=buf[k];
    else vector_atom[nlocal]=buf[0];
    return nvalues;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyPerAtom::pack_restart(int i, double *buf)
{
  buf[0] = nvalues+1;
  if (vectorStyle) for(int k=0;k<nvalues;k++) buf[k+1] = array_atom[i][k];
  else buf[0] = vector_atom[i];

  return (nvalues+1);
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyPerAtom::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[nlocal][k] = extra[nlocal][m++];
  else vector_atom[nlocal] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyPerAtom::maxsize_restart()
{
  return nvalues+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyPerAtom::size_restart(int nlocal)
{
  return nvalues+1;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int FixPropertyPerAtom::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
    int i,j;
    //we dont need to account for pbc here
    int m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      if (vectorStyle) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[j][k];
      else buf[m++] = vector_atom[j];
    }
    return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyPerAtom::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
      if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[i][k]=buf[m++];
      else vector_atom[i]=buf[m++];
  }

}

/* ---------------------------------------------------------------------- */

int FixPropertyPerAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (vectorStyle) for(int k=0;k<nvalues;k++) buf[m++] = array_atom[i][k];
    else buf[m++] = vector_atom[i];
  }
  return nvalues;
}

/* ---------------------------------------------------------------------- */

void FixPropertyPerAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (vectorStyle) for(int k=0;k<nvalues;k++) array_atom[j][k]+=buf[m++];
    else vector_atom[j]+=buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

