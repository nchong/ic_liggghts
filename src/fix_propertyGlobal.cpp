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
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "fix_propertyGlobal.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001
#define myAtof lmp->force->numeric 
#define myAtoi lmp->force->inumeric 

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::FixPropertyGlobal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    //Check args
    if (narg <6) error->all("Illegal fix property/global command, not enough arguments");

    //Read args
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    if (strcmp(arg[4],"scalar") == 0) svmStyle = FIXPROPERTYTYPE_GLOBAL_SCALAR;
    else if (strcmp(arg[4],"vector") == 0) svmStyle = FIXPROPERTYTYPE_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"peratomtype") == 0) svmStyle = FIXPROPERTYTYPE_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"matrix") == 0) svmStyle = FIXPROPERTYTYPE_GLOBAL_MATRIX;
    else if (strcmp(arg[4],"peratomtypepair") == 0) svmStyle = FIXPROPERTYTYPE_GLOBAL_MATRIX;
    else error->all("Unknown style for fix property/global. Valid styles are scalar or vector");

    int darg=0;
    if (svmStyle==FIXPROPERTYTYPE_GLOBAL_MATRIX) darg=1;

    //assign values to array
    nvalues = narg - 5 - darg;
    
    values = new double[nvalues];
    values_recomputed = new double[nvalues];
    for (int j=0;j<nvalues;j++) values[j] = myAtof(arg[5+darg+j]);

    if (svmStyle==FIXPROPERTYTYPE_GLOBAL_SCALAR) scalar_flag=1;
    else if (svmStyle==FIXPROPERTYTYPE_GLOBAL_VECTOR) {
        vector_flag=1;
        size_vector=nvalues;
    }
    else if (svmStyle==FIXPROPERTYTYPE_GLOBAL_MATRIX) {
        array_flag=1;
        size_array_cols=myAtoi(arg[5]);
        if (fmod(static_cast<double>(nvalues),size_array_cols)!=0.)
          error->all("Error in fix property/global: The number of default values must thus be a multiple of the nCols.");
        size_array_rows=static_cast<int>(static_cast<double>(nvalues)/size_array_cols);
    }

    extvector=0; 
    time_depend = 1; 

    //check if there is already a fix that tries to register a property with the same name
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyGlobal*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->all("Error in fix property/global. There is already a fix that registers a variable of the same name");

}

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::~FixPropertyGlobal()
{
  // delete locally stored arrays
  delete[] variablename;
  delete[] values;
}

/* ---------------------------------------------------------------------- */

bool FixPropertyGlobal::checkCorrectness(int desiredStyle, char* desiredVar,int desiredLen,int desiredRows)
{
    return ((desiredStyle==svmStyle)&&(!strcmp(desiredVar,variablename))&&(nvalues>=desiredLen)&&(svmStyle!=2||(svmStyle==2 && desiredRows>=size_array_rows)));

}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_scalar()
{
  return values[0];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_vector(int i)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    return values[i];
}

void FixPropertyGlobal::vector_modify(int i,double val)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    values_recomputed[i] = val;
}

double FixPropertyGlobal::compute_vector_modified(int i)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    return values_recomputed[i];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_array(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    int ind;
    if (i==0) ind=j;
    else ind =  ((i-1)*size_array_cols)+j;

    return values[ind];
}

void FixPropertyGlobal::array_modify(int i, int j,double val) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    int ind;
    if (i==0) ind=j;
    else ind =  ((i-1)*size_array_cols)+j;

    values_recomputed[ind] = val;
}

double FixPropertyGlobal::compute_array_modified(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    int ind;
    if (i==0) ind=j;
    else ind =  ((i-1)*size_array_cols)+j;

    return values_recomputed[ind];
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyGlobal::memory_usage()
{
  double bytes = nvalues * sizeof(double);
  return bytes;
}
