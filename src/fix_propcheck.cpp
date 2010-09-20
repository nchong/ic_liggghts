/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_propcheck.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

enum{PROP_FORCE,PROP_RADIUS};
enum{OP_GREATER,OP_SMALLER,OP_EQUAL};

/* ---------------------------------------------------------------------- */

FixPropCheck::FixPropCheck(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all("Illegal fix propcheck command");

  //scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 3;
  //global_freq = 1;
  //extscalar = 1;
  //extvector = 1;
  
  int iarg=3;

  nevery = atoi(arg[iarg++]);

  if(strcmp(arg[iarg++],"force")==0)
  {
     whichpropflag=PROP_FORCE;
  }
  //else if()...
  else error->all("Illegal fix propcheck command, 3rd argument");

  if(strcmp(arg[iarg++],">")==0)
  {
     opflag=OP_GREATER; //defined with enum above
  }  
  else if(strcmp(arg[iarg++],"<")==0)
  {
     opflag=OP_SMALLER;
  }  
  else if(strcmp(arg[iarg++],"=")==0)
  {
     opflag=OP_EQUAL;
  }  
  else error->all("Illegal fix propcheck command, 4th argument (check operator)");

  prop_threshold = atof(arg[iarg++]); // put an error message if of wrong type??

  if(strcmp(arg[iarg++],"radius")==0)
  {
     actionflag=PROP_RADIUS;
  }
  //else if()...
  else error->all("Illegal fix propcheck command, 6th argument");

  if(strcmp(arg[iarg++],"set")==0)
  {
     //maybe change lateron
  }
  //else if()...
  else error->all("Illegal fix propcheck command, 3rd argument");

  valtoset=atof(arg[iarg++]);

   // for lateron: possible use of set command
   /*char **setarg;
   setarg=new char*[4];
   for (int kk=0;kk<4;kk++) setarg[kk]=new char[30];
   setarg[0]="group";
   strcpy(setarg[1],);
   setarg[2]="neighlist/tri";
   strcpy(setarg[3],id);
   modify->add_fix(4,setarg);
   delete []fixarg;
   */
  
  fprintf(screen,"whichpropflag=%d, opflag=%d, actionflag=%d, prop_threshold=%f, valtoset=%f\n",whichpropflag,opflag,actionflag,prop_threshold ,valtoset);
  //error->all("george does not want to continue");

}

/* ---------------------------------------------------------------------- */

int FixPropCheck::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPropCheck::init()
{

}

/* ---------------------------------------------------------------------- */

void FixPropCheck::end_of_step()
{
  double *r = atom->radius;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double fmag;

  for (int i = 0; i < nlocal; i++) 
  {
    if (mask[i] & groupbit) {
      fmag = sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2]); // more efficient than squaring
      if (opflag == OP_GREATER){
        if (fmag > prop_threshold) { // can try doing it by variable.h?????
         r[i]=valtoset;
        }
      }
      else if (opflag == OP_SMALLER){
        if (fmag < prop_threshold) { // can try doing it by variable.h?????
         r[i]=valtoset;
        }
      }
      else if (opflag == OP_EQUAL) {
        if (fmag == prop_threshold) { // can try doing it by variable.h?????
         r[i]=valtoset;
        }
      }
      else error->all("What is the operator flag???");
    }
  }
}

/* ---------------------------------------------------------------------- */

