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

#include "compute_totalmass.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "stdlib.h"
#include "string.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTotalMass::ComputeTotalMass(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if ((narg != 3)&&(narg != 4)) error->all("Illegal compute totalmass command");

  filename=NULL;
  fp=NULL;

  if(narg==4)
  {
    filename=new char[strlen(arg[3]) + 1];
    strcpy(filename,arg[3]);
  }

  scalar_flag = 1;
  extscalar = 1;
  local_flag = 1;
  size_local_rows = 1;
  size_local_cols = 0;

  MPI_Comm_rank(world,&me);

  vector_local=new double[1];

  //erase old file
  if(me==0) fp = fopen(filename,"w");
  if(fp!=NULL) fclose(fp);

}

/* ---------------------------------------------------------------------- */

ComputeTotalMass::~ComputeTotalMass()
{
  delete []vector_local;
}

/* ---------------------------------------------------------------------- */

void ComputeTotalMass::init()
{
}

/* ---------------------------------------------------------------------- */

double ComputeTotalMass ::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  scalar = group->mass(igroup);
  if((me==0)&&(filename!=NULL)) fp = fopen(filename,"a");
  if(filename!=NULL)
  {
      fprintf(fp,"%d %f %f\n",update->ntimestep,static_cast<double>(update->ntimestep)*update->dt,scalar);
      fclose(fp);
  }
  return scalar;
}

void ComputeTotalMass ::compute_local()
{
  invoked_local = update->ntimestep;
  vector_local[0]=compute_scalar();
}
