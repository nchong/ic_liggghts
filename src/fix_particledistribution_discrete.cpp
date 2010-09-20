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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_sphere.h"
#include "fix_particledistribution_discrete.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  restart_global = 1;

  // random number generator, same for all procs

  if (narg < 7) error->all("Illegal fix particledistribution/discrete command, not enough arguments");
  seed=atoi(arg[3]);
  random = new RanPark(lmp,seed);
  ntemplates=atoi(arg[4]);
  if(ntemplates<1) error->all("Illegal fix particledistribution/discrete command, illegal number of ");

  templates=new FixTemplateSphere*[ntemplates];
  distweight=new double[ntemplates];
  cumweight=new double[ntemplates];
  parttogen=new int[ntemplates];
  distorder=new int[ntemplates];
  iarg=5;

  int itemp=0;

  //parse further args
  do{
    if(itemp==ntemplates) break;
    if(narg<iarg+1)error->all("Illegal fix particledistribution/discrete command, not enough arguments");
    int ifix=modify->find_fix(arg[iarg]);
    if(ifix<0) error->all("Illegal fix particledistribution/discrete command, invalid ID for fix particletemplate provided");
    if(strncmp(modify->fix[ifix]->style,"particletemplate/",16)!=0) error->all("Illegal fix particledistribution/discrete command, fix is not of type particletemplate");

    templates[itemp]=static_cast<FixTemplateSphere*>(modify->fix[ifix]);
    distweight[itemp]=atof(arg[iarg+1]);
    if (distweight[itemp] < 0) error->all("Illegal fix particledistribution/discrete command, invalid weight");
    itemp++;
    iarg += 2;
  }while (iarg < narg);

  //normalize distribution -
  double weightsum=0;
  for(int i=0;i<ntemplates;i++) weightsum+=distweight[i];
  if(fabs(weightsum-1.)>0.00001)error->warning("particledistribution/discrete: sum of distribution weights != 1, normalizing distribution");
  for(int i=0;i<ntemplates;i++) distweight[i]/=weightsum;

  if(comm->me==0&&screen){
      fprintf(screen,"particledistribution/discrete (id %s): distribution based on mass%%:\n",this->id);
      for(int i=0;i<ntemplates;i++)
        fprintf(screen,"    %s: d=%e (bounding sphere) mass%%=%f%%\n",templates[i]->id,2.*templates[i]->pti->r_bound,100.*distweight[i]);
  }

  //convert distribution from mass% to number%
  for(int i=0;i<ntemplates;i++) distweight[i]=distweight[i]/templates[i]->massexpect();
  weightsum=0;
  for(int i=0;i<ntemplates;i++) weightsum+=distweight[i];
  for(int i=0;i<ntemplates;i++) distweight[i]/=weightsum;

  if(comm->me==0&&screen){
      fprintf(screen,"particledistribution/discrete (id %s): distribution based on number%%:\n",this->id);
      for(int i=0;i<ntemplates;i++)
        fprintf(screen,"    %s: d=%e (bounding sphere) number%%=%f%%\n",templates[i]->id,2.*templates[i]->pti->r_bound,100.*distweight[i]);
  }

  cumweight[0]=distweight[0];
  for(int i=1;i<ntemplates;i++) cumweight[i]=distweight[i]+cumweight[i-1];

  volexpect=0.;massexpect=0.;
  for(int i=0;i<ntemplates;i++)
  {
      volexpect +=templates[i]->volexpect() *distweight[i];
      massexpect+=templates[i]->massexpect()*distweight[i];
  }

  //get min/maxtype
  maxtype=0;
  mintype=10000;
  for(int i=0;i<ntemplates;i++)
  {
    if(templates[i]->pti->atom_type>maxtype)
      maxtype=templates[i]->pti->atom_type;
    if(templates[i]->pti->atom_type<mintype)
      mintype=templates[i]->pti->atom_type;
  }

  maxnspheres=0;
  for(int i=0;i<ntemplates;i++)
    if(templates[i]->pti->nspheres>maxnspheres)
      maxnspheres=templates[i]->pti->nspheres;

  for(int i=0;i<ntemplates;i++) distorder[i]=i;

  bool swaped;
  int n=ntemplates;
  do
  {
      swaped=false;
      for(int i=0;i<ntemplates-1;i++)
      {
          if(templates[distorder[i]]->pti->r_bound<templates[distorder[i+1]]->pti->r_bound)
          {
            //swap
            int tmp=distorder[i];
            distorder[i]=distorder[i+1];
            distorder[i+1]=tmp;
            swaped=true;
          }
      }
      n--;
  }while(swaped&&n>0);

  pti=templates[distorder[0]]->pti;

}

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::~FixParticledistributionDiscrete()
{
    delete []templates;
    delete []distweight;
    delete []cumweight;
    delete []parttogen;
    delete []distorder;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::setmask()
{
  int mask = 0;
  return mask;
}

int FixParticledistributionDiscrete::random_init(int ntotal)
{
    ninsert=ntotal;
    ninserted=0;
    for(int i=0;i<ntemplates;i++) parttogen[i]=static_cast<int>(static_cast<double>(ninsert)*distweight[i]+random->uniform());

    int ntotal_new=0;
    for(int i=0;i<ntemplates;i++) ntotal_new+=parttogen[i];
    return ntotal_new;
}

/* ----------------------------------------------------------------------*/

void FixParticledistributionDiscrete::randomize()
{
    if(ntemplates==1){
         templates[0]->randomize();
         
         return; 
    }

    int chosen=0;
    int chosendist=distorder[chosen];
    int ntoinsert=parttogen[chosendist];
    while(ninserted>=ntoinsert&&chosen<ntemplates-1)
    {
        chosen++;
        chosendist=distorder[chosen];
        ntoinsert+=parttogen[chosendist];
    }

    templates[chosendist]->randomize();

    pti=templates[chosendist]->pti;

    ninserted++;

}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::vol_expect()
{
    return volexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::mass_expect()
{
    return massexpect;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::max_rad(int type)
{
    //get maxrad
    double maxrad=0.;
    for(int i=0;i<ntemplates;i++)
      if( templates[i]->pti->atom_type==type  && templates[i]->max_rad() > maxrad) maxrad = templates[i]->max_rad();

    return maxrad;
}
/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_type()
{
    return maxtype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::min_type()
{
    return mintype;
}

/* ----------------------------------------------------------------------*/

int FixParticledistributionDiscrete::max_nspheres()
{
    return maxnspheres;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);

  random->reset(seed);
}
