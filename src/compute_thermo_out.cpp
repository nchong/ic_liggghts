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

#include "compute_thermo_out.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "stdlib.h"
#include "string.h"
#include "modify.h"
#include "output.h"
#include "thermo.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeThermoOut::ComputeThermoOut(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal compute thermo/out command");

  fp_h=NULL;
  fp=NULL;

  filename=new char[strlen(arg[3]) + 1];
  strcpy(filename,arg[3]);

  filename_h=new char[strlen(arg[3]) + 8];
  sprintf(filename_h,"%s_header",filename);

  int iarg = 4;

  thermo_vals = 0;
  while (iarg<narg) thermo_vals |= (1 << atoi(arg[iarg++]));

  MPI_Comm_rank(world,&me);

  //erase old file
  if(me==0) fp = fopen(filename,"w");
  if(fp!=NULL) fclose(fp);
  if(me==0) fp_h = fopen(filename_h,"w");
  if(fp_h!=NULL) fclose(fp_h);

  if((me==0)&&(filename!=NULL))   fp   = fopen(filename,"a");
  if((me==0)&&(filename_h!=NULL)) fp_h = fopen(filename_h,"a");

  output->thermo->reg_th_cb(this);

  //some more init
  called_last=update->ntimestep-1;
  called_first=update->ntimestep-1;
  is_first_call=true;
}

/* ---------------------------------------------------------------------- */

ComputeThermoOut::~ComputeThermoOut()
{
    if(me==0 && fp!=NULL)fclose(fp);
    if(me==0 && fp_h!=NULL)fclose(fp_h);
    
    output->thermo->unreg_th_cb(this);
    delete []filename;
    delete []filename_h;
}

/* ---------------------------------------------------------------------- */

void ComputeThermoOut::init() {
  //print an extra newline at init
  
  if(me != 0) return;
  if(fp==NULL)   error->all("Compute thermo/out could not open file to output data");
  if(!is_first_call) fprintf(fp,"\n");
}

/* ---------------------------------------------------------------------- */

void ComputeThermoOut::callback(int ifield,bool isint,int intval,double doubleval,char *keyword)
{
  
  if(me != 0) return;

  if(!((1<<ifield) & thermo_vals)) return;

  if(fp==NULL) error->all("Compute thermo/out could not open file to output header");
  if(fp_h==NULL)   error->all("Compute thermo/out could not open file to output data");

  //check if first call
  if(is_first_call)
  {
      called_first=update->ntimestep;
      is_first_call=false;
  }

  //check if need to write header to file
  if(update->ntimestep == called_first) fprintf(fp_h,"%s,",keyword);

  //write numbers - write leading '\n' only if first write this ts and not first write
  if(update->ntimestep > called_last && update->ntimestep > called_first) fprintf(fp,"\n");

  if(isint) fprintf(fp,"%d\t",intval);
  else      fprintf(fp,"%f\t",doubleval);

  //store that has been called this ts
  called_last = update->ntimestep;
}

