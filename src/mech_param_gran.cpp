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
#include "lammps.h"
#include "atom.h"
#include "mpi.h"
#include "math.h"
#include "modify.h"
#include "mech_param_gran.h"
#include "error.h"
#include "fix_propertyGlobal.h"
#include "memory.h"
#include "fix_pour.h"
#include "fix_pour_dev.h"
#include "fix_wall_gran_hooke_history.h"
#include "fix_meshGran.h"

#define MUCH 1000000

using namespace LAMMPS_NS;

enum{MESHGRAN,XPLANE,YPLANE,ZPLANE,ZCYLINDER};

MechParamGran::MechParamGran(LAMMPS *lmp): Pointers(lmp)
{
    Yeff=NULL;
    veff=NULL;
    Geff=NULL;
    betaeff=NULL;
    cohEnergyDens=NULL;
    coeffRestLog=NULL;
    coeffFrict=NULL;
    charVel=0;
    arrays_active=false;
}

MechParamGran::~MechParamGran()
{
    destroy_arrays();
}

void MechParamGran::getMaterialParams(int charVelFlag, int cohesionflag)
{
  if(charVelFlag&&cohesionflag) error->warning("Cohesion model should only be used with hertzian contact laws.");

  //loop over all particles to check how many atom types are present
  min_type=1;
  max_type=1;

  for (int i=0;i<atom->nlocal;i++)
  {
      if (atom->type[i]<min_type) min_type=atom->type[i];
      if (atom->type[i]>max_type) max_type=atom->type[i];
  }

  //check all fixes of type pour
  for(int i=0;i<lmp->modify->nfix;i++)
  {
      
      if(strncmp(lmp->modify->fix[i]->style,"pour/dev",7)==0||strcmp(lmp->modify->fix[i]->style,"pour/multisphere")==0)
      {
          int tp_min=static_cast<FixPourDev*>(lmp->modify->fix[i])->ntype_min;
          int tp_max=static_cast<FixPourDev*>(lmp->modify->fix[i])->ntype_max;
          if(tp_min<min_type) min_type=tp_min;
          if(tp_max>max_type) max_type=tp_max;
      }
      else if(strncmp(lmp->modify->fix[i]->style,"pour",4)==0)
      {
          int tp=static_cast<FixPour*>(lmp->modify->fix[i])->ntype;
          if(tp<min_type) min_type=tp;
          if(tp>max_type) max_type=tp;
      }
      else if(strncmp(lmp->modify->fix[i]->style,"wall/gran",8)==0)
      {
          FixWallGranHookeHistory* fwg=static_cast<FixWallGranHookeHistory*>(lmp->modify->fix[i]);
          if(fwg->wallstyle==MESHGRAN)
          {
              for(int j=0;j<fwg->nFixMeshGran;j++)
              {
                int tp=fwg->FixMeshGranList[j]->atom_type_wall;
                if(tp<min_type) min_type=tp;
                if(tp>max_type) max_type=tp;
              }
          }
          else
          {
            int tp=fwg->atom_type_wall;
            if(tp<min_type) min_type=tp;
            if(tp>max_type) max_type=tp;
          }
      }
  }

  //Get min/max from other procs
  int min_type_all,max_type_all;
  MPI_Allreduce(&min_type,&min_type_all, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&max_type,&max_type_all, 1, MPI_INT, MPI_MAX, world);
  min_type=min_type_all;
  max_type=max_type_all;

  //error check
  if(min_type!=1) error->all("Atom types must start from 1 for granular simulations");
  if(max_type > atom->ntypes) error->all("Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");

  //Get pointer to the fixes that have the material properties
  
  Y1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0)]);
  v1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0)]);

  coeffRest1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type)]);
  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type)]);

  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type)]);
  if(charVelFlag)
    charVel1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0)]);

  if (arrays_active) destroy_arrays();
  create_arrays(max_type+1);

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          double Yi=Y1->compute_vector(i-1);
          double Yj=Y1->compute_vector(j-1);
          double vi=v1->compute_vector(i-1);
          double vj=v1->compute_vector(j-1);

          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);

          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);

          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
          //omitting veff here
      }
  }
  if(charVelFlag) charVel=charVel1->compute_scalar();

  //fprintf(lmp->screen,"Yeff[1][1]=%f Geff[1][1]=%f coeffRestLog[1][1]=%f coeffFrict[1][1]=%f betaeff[1][1]=%f\n",Yeff[1][1],Geff[1][1],coeffRestLog[1][1],coeffFrict[1][1],betaeff[1][1]);
}

void MechParamGran::create_arrays(int len)
{
    Yeff=memory->create_2d_double_array(len,len,"MechParamGran: Yeff");
    Geff=memory->create_2d_double_array(len,len,"MechParamGran: Geff");
    betaeff=memory->create_2d_double_array(len,len,"MechParamGran: Geff");
    veff=memory->create_2d_double_array(len,len,"MechParamGran: veff");
    cohEnergyDens=memory->create_2d_double_array(len,len,"MechParamGran: cohEnergyDens");
    coeffRestLog=memory->create_2d_double_array(len,len,"MechParamGran: coeffRest");
    coeffFrict=memory->create_2d_double_array(len,len,"MechParamGran: coeffFrict");
    arrays_active=true;
}

void MechParamGran::destroy_arrays()
{
    memory->destroy_2d_double_array(Yeff);
    memory->destroy_2d_double_array(Geff);
    memory->destroy_2d_double_array(betaeff);
    memory->destroy_2d_double_array(veff);
    memory->destroy_2d_double_array(cohEnergyDens);
    memory->destroy_2d_double_array(coeffRestLog);
    memory->destroy_2d_double_array(coeffFrict);
    arrays_active=false;
}
