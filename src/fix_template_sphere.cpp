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
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "domain.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixTemplateSphere::FixTemplateSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (domain->dimension != 3) error->all("Fix particletemplate/sphere is for 3D simulations only");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->all("Illegal fix particletemplate command, not enough arguments");
  seed=atoi(arg[3]);
  random = new RanPark(lmp,seed);
  PI = 4.0*atan(1.0);
  iarg = 4;

  if(strcmp(this->style,"particletemplate/sphere")==0)nspheres=1;
  else
  {
      nspheres=atoi(arg[iarg]);
      iarg++;
  }
  x_sphere=NULL;
  x_sphere_b=NULL;
  r_sphere=new double[nspheres];

  //set default values
  atom_type=1;
  density1=1.;
  density2=0.;
  density_randstyle=RAN_STYLE_CONSTANT_FTS;
  r_sphere[0]=1.;
  radius2=1.;
  radius_randstyle=RAN_STYLE_CONSTANT_FTS;

  pti=new ParticleToInsert(lmp,nspheres);

  //parse further args
  while (iarg < narg) {
    if (strcmp(arg[iarg],"atom_type") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
      atom_type=atoi(arg[iarg+1]);
      if (atom_type < 1) error->all("Illegal fix particletemplate command, invalid atom type (must be >=1)");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
      if (strcmp(arg[iarg+1],"constant")==0) density_randstyle=RAN_STYLE_CONSTANT_FTS;
      else if (strcmp(arg[iarg+1],"uniform")==0) density_randstyle=RAN_STYLE_UNIFORM_FTS;
      else if (strcmp(arg[iarg+1],"gaussian")==0) density_randstyle=RAN_STYLE_GAUSSIAN_FTS;

      density1=atof(arg[iarg+2]);
      if(density_randstyle!=RAN_STYLE_CONSTANT_FTS)
      {
          if (iarg+3 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
          density2=atof(arg[iarg+3]);
          iarg++;
      }
      if (density1 < 0 || density2 < 0) error->all("Illegal fix particletemplate command, invalid density");
      iarg += 3;
    }
    else if (strcmp(arg[iarg],"radius") == 0) {
      if(strcmp(this->style,"particletemplate/sphere")!=0) error->all("Illegal fix particletemplate command, keyword radius only valid for particletemplate/sphere");
      if (iarg+3 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
      if (strcmp(arg[iarg+1],"constant")==0) radius_randstyle=RAN_STYLE_CONSTANT_FTS;
      else if (strcmp(arg[iarg+1],"uniform")==0) radius_randstyle=RAN_STYLE_UNIFORM_FTS;
      else if (strcmp(arg[iarg+1],"gaussian")==0) radius_randstyle=RAN_STYLE_GAUSSIAN_FTS;

      r_sphere[0]=atof(arg[iarg+2]);
      if(r_sphere[0]<=0.) error->all("Illegal fix particletemplate/sphere command, radius must be >0");
      if(radius_randstyle!=RAN_STYLE_CONSTANT_FTS)
      {
          if (iarg+3 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
          radius2=atof(arg[iarg+3]);
          if(radius_randstyle==RAN_STYLE_UNIFORM_FTS) radius2=radius2/r_sphere[0];
          iarg++;
      }
      if (r_sphere[0] < 0) error->all("Illegal fix particletemplate command, invalid radius");
      iarg += 3;
    }
    else if (strcmp(arg[iarg],"spheres")==0||strcmp(arg[iarg],"ntry")==0) break;
    else error->all("Illegal fix particletemplate command, unrecognized keyword");
  }

  if(strcmp(this->style,"particletemplate/sphere")!=0)return;

  //set mass and volume
  volume=pow((2.*r_sphere[0]),3.)*PI/6.;
  mass=density1*volume;

  //copy the values to public
  pti->nspheres=nspheres;
  pti->density_ins=density1;
  pti->volume_ins=volume;
  pti->mass_ins=mass;
  pti->radius_ins[0]=r_sphere[0];
  pti->r_bound=r_sphere[0];
  for (int i=0;i<3;i++) pti->x_ins[0][i]=0.;

  pti->atom_type=atom_type;
}

/* ---------------------------------------------------------------------- */

FixTemplateSphere::~FixTemplateSphere()
{
    delete []r_sphere;
    delete random;
    delete pti;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::randomize()
{
   randomize_r();
   randomize_dens();
   
}

/* ----------------------------------------------------------------------*/

inline void FixTemplateSphere::randomize_r()
{
   if (radius_randstyle==RAN_STYLE_CONSTANT_FTS) return;

   //use distribution based on mass-%
   if (radius_randstyle==RAN_STYLE_UNIFORM_FTS)
   {
       double rn = random->uniform();
       rn = rn*rn*rn;
       pti->radius_ins[0]=r_sphere[0]+(radius2*r_sphere[0]-r_sphere[0])*rn;
   }
   else if (radius_randstyle==RAN_STYLE_GAUSSIAN_FTS)
   {
       error->all("Random style gaussian not available for radius");//how to do transformation to mass %??
       double rn = random->gaussian();
       pti->radius_ins[0]=r_sphere[0]+(radius2*r_sphere[0]-r_sphere[0])*rn;
       if(pti->radius_ins[0]<0)
       {
           pti->radius_ins[0]=r_sphere[0]/10.;
           error->warning("Insertion radius value generated is <0, assuming some value");
       }
   }
   pti->r_bound=pti->radius_ins[0];
}

double FixTemplateSphere::max_rad()
{
    double maxrad=0.;
    for(int j=0;j<nspheres;j++)
          if(r_sphere[j]>maxrad) maxrad=r_sphere[j];

    if (radius_randstyle==RAN_STYLE_CONSTANT_FTS) return maxrad;
    if (radius_randstyle==RAN_STYLE_UNIFORM_FTS) return radius2*maxrad;
    if (radius_randstyle==RAN_STYLE_GAUSSIAN_FTS) error->all("Random style gaussian not available for radius");
    return 0.;

}

/* ----------------------------------------------------------------------*/

inline void FixTemplateSphere::randomize_dens()
{
   if (density_randstyle==RAN_STYLE_CONSTANT_FTS) return;
   else if (density_randstyle==RAN_STYLE_UNIFORM_FTS)
   {
       pti->density_ins=density1+(density2-density1)*(random->uniform());
   }
   else if (density_randstyle==RAN_STYLE_GAUSSIAN_FTS)
   {
       pti->density_ins=density1+density2*(random->gaussian());
       if(pti->density_ins<0)
       {
           pti->density_ins=density1/10.;
           error->warning("density value generated is <0, assuming some value");
       }
   }

   pti->volume_ins=pow((2.*pti->radius_ins[0]),3.)*PI/6.;
   pti->mass_ins=pti->density_ins*pti->volume_ins;
}

/* ----------------------------------------------------------------------*/
//volume expectancy for different random styles
//uniform:  param 1 = low  param 2 = high
//gaussian: param 1 = mu   param 2 = sigma

double FixTemplateSphere::volexpect()
{
    double param_1=r_sphere[0];
    double param_2=radius2;
    if (radius_randstyle==RAN_STYLE_UNIFORM_FTS) param_2=r_sphere[0]*radius2;

    if(radius_randstyle==RAN_STYLE_CONSTANT_FTS)
      return volume;

    if(radius_randstyle==RAN_STYLE_UNIFORM_FTS)
      return (PI/3.0 * (pow(param_2,4.)-pow(param_1,4.))/(param_2-param_1));

    if(radius_randstyle==RAN_STYLE_GAUSSIAN_FTS)
      return (4.*PI/3.*param_1*(param_1*param_1+3.*param_2*param_2));
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::massexpect()
{
    if      (density_randstyle==RAN_STYLE_CONSTANT_FTS)
      return density1*volexpect();

    else if (density_randstyle==RAN_STYLE_UNIFORM_FTS)
      return 0.5*(density1+density2)*volexpect();

    else if (density_randstyle==RAN_STYLE_GAUSSIAN_FTS)
      return density1*volexpect();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTemplateSphere::write_restart(FILE *fp)
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

void FixTemplateSphere::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);

  random->reset(seed);
}
