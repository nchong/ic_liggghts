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
#include "comm.h"
#include "myvector.h"
#include "fix_cfd_coupling.h"
#include "fix_propertyGlobal.h"
#include "fix_propertyPerAtom.h"
#include "fix_cfd_coupling.h"
#include "cfd_regionmodel.h"
#include "style_cfd_datacoupling.h"
#include "style_cfd_regionmodel.h"

using namespace LAMMPS_NS;

#define MAXLENGTH 30

/* ---------------------------------------------------------------------- */

FixCfdCoupling::FixCfdCoupling(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  for(int ii=0;ii<modify->nfix;ii++)
  {
      if(strcmp(modify->fix[ii]->style,style) == 0) error->all("There must not be more than one fix of type couple/cfd");
  }

  if (narg < 5) error->all("Illegal fix couple/cfd/file command");

  int iarg = 3;
  couple_nevery = atoi(arg[iarg++]);
  if(couple_nevery<=0)error->all("Fix couple/cfd/file: nevery must be >0");

  nevery = 1;

  if (0) return;
  #define CFD_DATACOUPLING_CLASS
  #define CfdDataCouplingStyle(key,Class) \
  else if (strcmp(arg[iarg],#key) == 0) dc = new Class(lmp,iarg+1,narg,arg,this);
  #include "style_cfd_datacoupling.h"
  #undef CFD_DATACOUPLING_CLASS
  else error->all("Unknown cfd data coupling style");

  iarg = dc->get_iarg();

  rm = NULL;

  if(iarg < narg)
  {
      if (0) return;
      #define CFD_REGIONMODEL_CLASS
      #define CfdRegionStyle(key,Class) \
      else if (strcmp(arg[iarg],#key) == 0) rm = new Class(lmp,iarg+1,narg,arg,this);
      #include "style_cfd_regionmodel.h"
      #undef CFD_REGIONMODEL_CLASS
      else error->all("Unknown cfd regionmodel style");
  }

  if(rm) iarg = rm->get_iarg();

  extvector = 1; //extensive
  create_attribute = 1;

  peratom_flag = 1;
  peratom_freq=1;

  array_flag = 1;
  array_atom = NULL;

  atom->add_callback(0);

  npull = 0;
  npush = 0;
  nvalues_max = 0;

  pullnames = NULL;
  pulltypes = NULL;
  pushnames = NULL;
  pushtypes = NULL;
  grow_();

  ts_create = update->ntimestep;
}

void FixCfdCoupling::grow_()
{
      nvalues_max+=10;
      pullnames = memory->grow_2d_char_array(pullnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:valnames");
      pulltypes = memory->grow_2d_char_array(pulltypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:valtypes");
      pushnames = memory->grow_2d_char_array(pushnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushnames");
      pushtypes = memory->grow_2d_char_array(pushtypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushtypes");
}

FixCfdCoupling::~FixCfdCoupling()
{
	atom->delete_callback(id,0);
	memory->destroy_2d_double_array(array_atom);
	memory->destroy_2d_char_array(pullnames);
	memory->destroy_2d_char_array(pulltypes);
	memory->destroy_2d_char_array(pushnames);
	memory->destroy_2d_char_array(pushtypes);
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

void FixCfdCoupling::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if(rm) rm->init();

  init_submodel();
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

void FixCfdCoupling::add_pull_property(char *name,char *type)
{
    if(strlen(name) >= MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");
    if(npull >= nvalues_max) grow_();

    for(int i = 0; i < npull; i++)
    {
        if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type) == 0) return;
        if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type)) error->all("Properties added via FixCfdCoupling::add_pull_property are inconsistent");
    }

    strcpy(pullnames[npull],name);
    strcpy(pulltypes[npull],type);
    npull++;
}

/* ---------------------------------------------------------------------- */
void FixCfdCoupling::add_push_property(char *name,char *type)
{
    if(strlen(name) >= MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");
    if(npush >= nvalues_max) grow_();

    for(int i = 0; i < npush; i++)
    {
        if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type) == 0) return;
        if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type)) error->all("Properties added via FixCfdCoupling::add_push_property are inconsistent");
    }

    strcpy(pushnames[npush],name);
    strcpy(pushtypes[npush],type);
    npush++;
}

/* ---------------------------------------------------------------------- */

void* FixCfdCoupling::find_pull_property(char *name,char *type,int &len1,int &len2)
{
    return find_property(0,name,type,len1,len2);
}
void* FixCfdCoupling::find_push_property(char *name,char *type,int &len1,int &len2)
{
    return find_property(1,name,type,len1,len2);
}

void* FixCfdCoupling::find_property(int push,char *name,char *type,int &len1,int &len2)
{
    fprintf(screen,"find_property called for %s, var %s, type %s\n",push==1?"push":"pull",name,type);

    void *ptr = NULL;
    int flag = 0;

    //check existence
    if(push)
    {
        for(int i = 0; i < npush; i++)
            if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type) == 0) flag = 1;
    }
    else
    {
        for(int i = 0; i < npull; i++)
            if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type) == 0) flag = 1;
    }

    if(!flag) error->all("Inconsistency in FixCfdCoupling::find_property");

    if(atom->extract(name)) return atom->extract(name);

    int ifix1 = -1, ifix2 = -1;

    if(strcmp(type,"scalar") == 0 || strcmp(type,"vector") == 0)
       ifix1 = modify->find_fix_property(name,"property/peratom",type,0,0,false);
    else if(strcmp(type,"globalscalar") == 0)
       ifix2 = modify->find_fix_property(name,"property/global","scalar",0,0,false);
    else if(strcmp(type,"globalvector") == 0)
       ifix2 = modify->find_fix_property(name,"property/global","vector",0,0,false);
    else if(strcmp(type,"globalmatrix") == 0)
       ifix2 = modify->find_fix_property(name,"property/global","matrix",0,0,false);

    if(ifix1 > -1 && strcmp(type,"scalar") == 0) ptr = (void*) static_cast<FixPropertyPerAtom*>(modify->fix[ifix1])->vector_atom;
    if(ifix1 > -1 && strcmp(type,"vector") == 0) ptr = (void*) static_cast<FixPropertyPerAtom*>(modify->fix[ifix1])->array_atom;

    if(ifix2 > -1 && strcmp(type,"globalvector") == 0)
    {
        ptr = (void*) static_cast<FixPropertyGlobal*>(modify->fix[ifix2])->values;
        len1 = static_cast<FixPropertyGlobal*>(modify->fix[ifix2])->nvalues;
    }

    if(ifix2 > -1 && strcmp(type,"globalarray") == 0)
    {
        ptr = (void*) static_cast<FixPropertyGlobal*>(modify->fix[ifix2])->array;
        len1  = static_cast<FixPropertyGlobal*>(modify->fix[ifix2])->size_array_rows;
        len2  = static_cast<FixPropertyGlobal*>(modify->fix[ifix2])->size_array_cols;
    }
    return ptr;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::end_of_step()
{
    
    if(!dc->liggghts_is_active) return;

    int ts = update->ntimestep;

    if((ts+1) % couple_nevery || ts_create == ts+1) couple_this = 0;
    else couple_this = 1;

    if(ts % couple_nevery || ts_create == ts) return;

    if(screen && comm->me == 0) fprintf(screen,"CFD Coupling established at step %d\n",ts);

    void *dummy = NULL;

    if(ts % couple_nevery == 0)
    {
      //call region model
      if(rm) rm->rm_update();

      //write to file
      for(int i = 0; i < npush; i++)
      {
           
           dc->push(pushnames[i],pushtypes[i],dummy);
      }

      //read from files
      for(int i = 0; i < npull; i++)
      {
         
         dc->pull(pullnames[i],pulltypes[i],dummy);
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
