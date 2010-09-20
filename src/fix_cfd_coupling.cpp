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
#include "fix_cfd_coupling.h"
#include "style_datacoupling.h"

using namespace LAMMPS_NS;

#define MAXLENGTH 30

/* ---------------------------------------------------------------------- */

FixCfdCoupling::FixCfdCoupling(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  for(int ii=0;ii<modify->nfix;ii++)
  {
      if(strcmp(modify->fix[ii]->style,style)==0) error->all("There must not be more than one fix of type couple/cfd");
  }

  if (narg < 5) error->all("Illegal fix couple/cfd/file command");

  int iarg = 3;
  couple_nevery = atoi(arg[iarg++]);
  if(couple_nevery<=0)error->all("Fix couple/cfd/file: nevery must be >0");

  nevery = 1;

  if (0) return;
  #define COUPLING_CLASS
  #define DataCouplingStyle(key,Class) \
  else if (strcmp(arg[iarg],#key) == 0) dc = new Class(lmp,narg-iarg-1,&arg[iarg+1],this);
  #include "style_datacoupling.h"
  #undef COUPLING_CLASS
  else error->all("Unknown data coupling style");

  extvector = 1; //extensive
  create_attribute = 1;

  peratom_flag = 1;
  peratom_freq=1;

  array_flag = 1;
  array_atom = NULL;

  atom->add_callback(0);

  nvalues = 0;
  npush = 0;
  nvalues_max = 0;

  valnames = NULL;
  valtypes = NULL;
  nvals_each = NULL;
  pushnames = NULL;
  pushtypes = NULL;
  grow_();
}

void FixCfdCoupling::grow_()
{
      nvalues_max+=10;
      valnames = memory->grow_2d_char_array(valnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:valnames");
      valtypes = memory->grow_2d_char_array(valtypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:valtypes");
      nvals_each = (int*)memory->srealloc(nvals_each,nvalues_max*sizeof(int),"FixCfdCoupling:nvals_each");
      pushnames = memory->grow_2d_char_array(pushnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushnames");
      pushtypes = memory->grow_2d_char_array(pushtypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushtypes");
}

FixCfdCoupling::~FixCfdCoupling()
{
	atom->delete_callback(id,0);
	memory->destroy_2d_double_array(array_atom);
	memory->destroy_2d_char_array(valnames);
	memory->destroy_2d_char_array(valtypes);
}

/* ---------------------------------------------------------------------- */

int FixCfdCoupling::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::set_arrays(int i)
{
    for (int m = 0; m < nvalues; m++)
		array_atom[i][m] = 0.;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::setup(int vflag)
{

  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }

  if(update->ntimestep == 0) end_of_step();
}

/* ---------------------------------------------------------------------- */
void FixCfdCoupling::push(char *name,char *type,void *&ptr)
{
    return dc->push(name,type,ptr);
}

void FixCfdCoupling::pull(char *name,char *type,void *&ptr)
{
    return dc->pull(name,type,ptr);
}

/* ---------------------------------------------------------------------- */

int FixCfdCoupling::add_pull_property(char *name,char *type)
{
    if(strlen(name)>=MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");

    if(nvalues>=nvalues_max) grow_();

    strcpy(valnames[nvalues],name);
    strcpy(valtypes[nvalues],type);

    if (strcmp(type,"vector")==0) nvals_each[nvalues] = 3;
    else if (strcmp(type,"scalar")==0) nvals_each[nvalues] = 1;
    else error->all("ix couple/cfd: Unregognized data type");
    nvalues++;

    //return index in array_atom where it is located
    int index = 0;
    for(int i=0;i<nvalues-1;i++)
      index += nvals_each[i];

    return index;
}

/* ---------------------------------------------------------------------- */

void* FixCfdCoupling::find_pull_property(char *name,char *type)
{
    int iprop = -1;
    for(int i=0;i<nvalues;i++)
        if(strcmp(name,valnames[i]) ==0 && strcmp(type,valtypes[i])==0) iprop = i;

    if(iprop == -1) error->all("Fix couple/cfd could not locate pull property");

    int index = 0;
    for(int i=0;i<iprop;i++)
      index += nvals_each[i];

    return (void*)(&array_atom[index]);
}

/* ---------------------------------------------------------------------- */
void FixCfdCoupling::add_push_property(char *name,char *type)
{
    if(strlen(name)>=MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");
    if(npush>=nvalues_max) grow_();

    strcpy(pushnames[npush],name);
    strcpy(pushtypes[npush],type);
    npush++;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::end_of_step()
{
    
    if(!dc->liggghts_is_active) return;

    int ts = update->ntimestep;

    if(ts % couple_nevery) return;

    void *dummy = NULL;

    if(ts % couple_nevery == 0)
    {

      //write to file  
      for(int i = 0;i<npush;i++)
      {
           
           dc->push(pushnames[i],pushtypes[i],dummy);
      }

      //read from files
      for(int i = 0;i<nvalues;i+=nvals_each[i])
      {
         
         dc->pull(valnames[i],valtypes[i],dummy);
      }

    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return stored values
------------------------------------------------------------------------- */

double FixCfdCoupling::compute_array(int i,int j)
{
    
    return array_atom[i][j];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixCfdCoupling::memory_usage()
{
  double bytes;
  bytes = atom->nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCfdCoupling::grow_arrays(int nmax)
{
    array_atom = memory->grow_2d_double_array(array_atom,nmax,nvalues,"fix_cfd_coupling:array_atom");
}

void FixCfdCoupling::zero_arrays()
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
	set_arrays(i);
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCfdCoupling::copy_arrays(int i, int j)
{
  for (int m = 0; m < nvalues; m++)
    array_atom[j][m] = array_atom[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCfdCoupling::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = array_atom[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCfdCoupling::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) array_atom[nlocal][m] = buf[m];
  return nvalues;
}

