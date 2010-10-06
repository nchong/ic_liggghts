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
#include "fix_heat_gran.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "math_extra.h"
#include "fix_propertyGlobal.h"
#include "fix_propertyPerAtom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixHeatGran::FixHeatGran(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all("Fix heat/gran needs per particle radius and mass");

  //check if a fix of this style already exists
  for (int i=0;i<modify->nfix;i++)
      if (strcmp(modify->fix[i]->style,style) == 0) error->all("A fix of type heat/gran is already registered. Can not have more than one");

  if (narg < 4) error->all("Illegal fix heat/gran command, not enough arguments");
  T0 = atof(arg[3]);

  fppat = static_cast<FixPropertyPerAtom*>(NULL);
  fppahf = static_cast<FixPropertyPerAtom*>(NULL);
  fppahs = static_cast<FixPropertyPerAtom*>(NULL);

  peratom_flag=1;              
  size_peratom_cols=0;         
  peratom_freq=1;
  time_depend=1;

  scalar_flag=1; 
  global_freq=1; 

}

/* ---------------------------------------------------------------------- */

FixHeatGran::~FixHeatGran()
{
    
    //unregister Temp as property/peratom
    if (fppat!=NULL) modify->delete_fix("Temp");
    //unregister heatFlux as property/peratom
    if (fppahf!=NULL) modify->delete_fix("heatFlux");
    //unregister heatSource as property/peratom
    if (fppahs!=NULL) modify->delete_fix("heatSource");

}

/* ---------------------------------------------------------------------- */

int FixHeatGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  return mask;
}
/* ---------------------------------------------------------------------- */
void FixHeatGran::updatePtrs()
{
  Temp=fppat->vector_atom;
  vector_atom=Temp; 

  heatFlux=fppahf->vector_atom;

  heatSource=fppahs->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::init()
{
  if (!atom->radius_flag) error->all("Please use a granular pair style for fix heat/gran");
  if (!atom->rmass_flag) error->all("Please use a granular pair style for fix heat/gran");

  char **fixarg;
  fixarg=new char*[9];
  for (int kk=0;kk<9;kk++) fixarg[kk]=new char[30];

  if (fppat==NULL) {
  //register Temp as property/peratom
    fixarg[0]="Temp";
    fixarg[1]="all";
    fixarg[2]="property/peratom";
    fixarg[3]="Temp";
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    sprintf(fixarg[8],"%f",T0);
    modify->add_fix(9,fixarg);
    fppat=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("Temp","property/peratom","scalar",0,0)]);
  }

  if (fppahf==NULL){
    //register heatFlux as property/peratom
    fixarg[0]="heatFlux";
    fixarg[1]="all";
    fixarg[2]="property/peratom";
    fixarg[3]="heatFlux";
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="no";    
    fixarg[7]="yes";    
    fixarg[8]="0.";     
    modify->add_fix(9,fixarg);
    fppahf=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("heatFlux","property/peratom","scalar",0,0)]);
  }

  if (fppahs==NULL){
    //register heatSource as property/peratom
    fixarg[0]="heatSource";
    fixarg[1]="all";
    fixarg[2]="property/peratom";
    fixarg[3]="heatSource";
    fixarg[4]="scalar"; 
    fixarg[5]="yes";    
    fixarg[6]="yes";    
    fixarg[7]="no";    
    fixarg[8]="0.";     
    modify->add_fix(9,fixarg);
    fppahs=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("heatSource","property/peratom","scalar",0,0)]);
  }
  delete []fixarg;

  int nAtomTypes=atom->ntypes;

  fpgco=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("thermalConductivity","property/global","peratomtype",atom->ntypes,0)]);
  fpgca=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("thermalCapacity","property/global","peratomtype",atom->ntypes,0)]);

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_force(int vflag)
{
  double hc,contactArea;
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,tcoi,tcoj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0) error->warning("Fix heat/gran may implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0) error->warning("Fix heat/gran may implementation may not be valid for pair style hybrid/overlay");

  inum = force->pair->list->inum;
  ilist = force->pair->list->ilist;
  numneigh = force->pair->list->numneigh;
  firstneigh = force->pair->list->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();

  //reset heat flux
  for (int i = 0; i < atom->nlocal; i++)
  {
       if (mask[i] & groupbit){
           heatFlux[i]=0.;
       }
  }

  fppat->do_forward_comm();
  fppahs->do_forward_comm();

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq < radsum*radsum) {  //contact
         r = sqrt(rsq);
         contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r); //contact area of the two spheres
         tcoi=fpgco->compute_vector(type[i]-1); //types start at 1, array at 0
         tcoj=fpgco->compute_vector(type[j]-1);

         if ((fabs(tcoi)<1e-7)||(fabs(tcoj)<1e-7)) hc=0.;
         else hc=4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

         heatFlux[i]+=(Temp[j]-Temp[i])*hc+heatSource[i];
         if (newton_pair||j<nlocal) heatFlux[j]+=(Temp[i]-Temp[j])*hc+heatSource[j];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixHeatGran::final_integrate()
{
    double dt=update->dt;
    int nlocal = atom->nlocal;
    double *rmass = atom->rmass;
    int *type = atom->type;
    int *mask = atom->mask;
    double tcai;

    updatePtrs();

    for (int i = 0; i < nlocal; i++)
    {
       if (mask[i] & groupbit){
          tcai=fpgca->compute_vector(type[i]-1);
          if(fabs(tcai)>1.e-3) Temp[i]+=heatFlux[i]*dt/(rmass[i]*tcai);
       }
    }
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  error->all("Can not use fix heat/gran together with rRESPA");
}

/* ---------------------------------------------------------------------- */
double FixHeatGran::compute_scalar()
{

    double *rmass = atom->rmass;
    int *type = atom->type;
    double tcai;

    updatePtrs();

    double heat_energy=0.;
    for (int i = 0; i < atom->nlocal; i++)
    {
       tcai=fpgca->compute_vector(type[i]-1);

       heat_energy+=tcai*rmass[i]*Temp[i];
    }
    double heat_energy_all=0.;

    MPI_Allreduce(&heat_energy,&heat_energy_all,1,MPI_DOUBLE, MPI_SUM, world);
    return heat_energy_all;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::reset_dt()
{
    error->all("Adaptive time-stepping not implemented for fix heat/gran");
}
