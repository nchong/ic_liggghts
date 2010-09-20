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
#include "fix_meshGran_analyze.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "triSpherePenetration.h"
#include "stl_tri.h"
#include "mpi.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001

#define DEBUGMODE_MESHGRANANALYZE false
#define DEBUG_OUTP_MESHGRAN logfile

/* ---------------------------------------------------------------------- */

FixMeshGranAnalyze::FixMeshGranAnalyze(LAMMPS *lmp, int narg, char **arg) :
  FixMeshGran(lmp, narg, arg)
{
   tmp=new double[3];
   tmp2=new double[3];
   analyseStress=true;

   vector_flag = 1;
   size_vector = 6;
   global_freq = 1;
   extvector = 1;

}

FixMeshGranAnalyze::~FixMeshGranAnalyze()
{
   delete []tmp;
   delete []tmp2;
}

/* ---------------------------------------------------------------------- */

int FixMeshGranAnalyze::setmask()
{
  int mask = 0;
  mask |=PRE_FORCE;
  mask |=FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::pre_force(int vflag)
{
    force_total[0]=0.;  force_total[1]=0.;  force_total[2]=0.;
    torque_total[0]=0.; torque_total[1]=0.; torque_total[2]=0.;

    for(int i=0;i<nTri;i++)
        for(int j=0;j<3;j++)
            STLdata->f_tri[i][j]=0.;

}

/* ---------------------------------------------------------------------- */
//   called during wall force calc
/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::add_particle_contribution(double *frc,double* contactPoint,int iTri)
{
    
    vectorAdd3D(STLdata->f_tri[iTri],frc,STLdata->f_tri[iTri]);

    vectorAdd3D(force_total,frc,force_total);
    vectorSubtract3D(contactPoint,p_ref,tmp);
    vectorCross3D(tmp,frc,tmp2); //tmp2 is torque contrib
    vectorAdd3D(torque_total,tmp2,torque_total);
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::final_integrate()
{
    calc_total_force();
}

/* ---------------------------------------------------------------------- */

void FixMeshGranAnalyze::calc_total_force()
{
    //total force on tri
    double *force_total_all=new double[3];
    double *torque_total_all=new double[3];

    MPI_Allreduce(force_total,force_total_all,3, MPI_DOUBLE, MPI_SUM,world);
    MPI_Allreduce(torque_total,torque_total_all,3, MPI_DOUBLE, MPI_SUM,world);

    for(int j=0;j<3;j++)
    {
        force_total[j]=force_total_all[j];
        torque_total[j]=torque_total_all[j];
    }

    delete []force_total_all;
    delete []torque_total_all;

    //forces on tri

    MPI_Allreduce(&(STLdata->f_tri[0][0]),&(STLdata->fn_fshear[0][0]),3*(STLdata->nTri), MPI_DOUBLE, MPI_SUM,world);

    //switch fn_fshear and f_tri
    double **helper;
    helper=STLdata->f_tri;
    STLdata->f_tri=STLdata->fn_fshear;
    STLdata->fn_fshear=helper;

    double *temp=new double[3];
    double p,s;
    for(int i=0;i<nTri;i++)
    {
        //pressure
        STLdata->fn_fshear[i][0]=vectorDot3D(STLdata->f_tri[i],STLdata->facenormal[i]);
        vectorScalarMult3D(STLdata->facenormal[i],STLdata->fn_fshear[i][0],temp);
        vectorSubtract3D(STLdata->f_tri[i],temp,temp);
        //shear force
        STLdata->fn_fshear[i][1]=vectorMag3D(temp);
    }
    delete []temp;
}

/* ----------------------------------------------------------------------
   return force/torque on body
------------------------------------------------------------------------- */

double FixMeshGranAnalyze::compute_vector(int n)
{
  //return force/torque
  if(n<3) return force_total[n];
  else    return torque_total[n-3];
}
