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
#include "stdio.h"
#include "fix_shear_history.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define DELTA_MAXTOUCH 15

/* ---------------------------------------------------------------------- */

FixShearHistory::FixShearHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  restart_peratom = 1;
  create_attribute = 1;

  maxtouch=DELTA_MAXTOUCH;

  // perform initial allocation of atom-based arrays
  // register with atom class

  npartner = NULL;
  partner = NULL;
  shearpartner = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

FixShearHistory::~FixShearHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->sfree(npartner);
  memory->destroy_2d_int_array(partner);
  memory->destroy_3d_double_array(shearpartner);
}

/* ---------------------------------------------------------------------- */

int FixShearHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShearHistory::init()
{
  if (atom->tag_enable == 0)
    error->all("Pair style granular with history requires atoms have IDs");
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be exchanged with atoms
------------------------------------------------------------------------- */

void FixShearHistory::pre_exchange()
{
  int i,j,ii,jj,m,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  // zero npartners for all current atoms

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) npartner[i] = 0;

  // copy shear info from neighbor list atoms to atom arrays

  int *tag = atom->tag;

  NeighList *list = pair->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listgranhistory->firstneigh;
  firstshear = list->listgranhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allshear = firstshear[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
	shear = &allshear[3*jj];
	j = jlist[jj];
	  if (npartner[i] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 
	  
	  m = npartner[i];
	  partner[i][m] = tag[j];
	  shearpartner[i][m][0] = shear[0];
	  shearpartner[i][m][1] = shear[1];
	  shearpartner[i][m][2] = shear[2];
	  
	npartner[i]++;
	if (j < nlocal) {
        if (npartner[j] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 
	    
	    m = npartner[j];
	    partner[j][m] = tag[i];
	    shearpartner[j][m][0] = -shear[0];
	    shearpartner[j][m][1] = -shear[1];
	    shearpartner[j][m][2] = -shear[2];
	    
	  npartner[j]++;
	}
      }
    }
  }

  int maxtouch_all;
  MPI_Allreduce(&maxtouch,&maxtouch_all,1,MPI_INT,MPI_MAX,world);
  while (maxtouch<maxtouch_all) grow_arrays_maxtouch(atom->nmax);

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixShearHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*maxtouch * sizeof(int);  
  bytes += nmax*maxtouch *3 * sizeof(double); 
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::grow_arrays(int nmax)
{
  npartner = (int *) memory->srealloc(npartner,nmax*sizeof(int),
				      "shear_history:npartner");
  partner = memory->grow_2d_int_array(partner,nmax,maxtouch,
				      "shear_history:partner");
  //fprintf(screen,"nmax=%d, maxtouch=%d\n",nmax,maxtouch);
  shearpartner =
    memory->grow_3d_double_array(shearpartner,nmax,maxtouch,3,
				 "shear_history:shearpartner");
}

/* ----------------------------------------------------------------------
   grow local atom-based arrays in case maxtouch is too small 
------------------------------------------------------------------------- */

void FixShearHistory::grow_arrays_maxtouch(int nmax)
{
  if(comm->me==0)
  {
      if(screen) fprintf(screen,"INFO: more than %d touching neighbor atoms found, growing shear history...",maxtouch);
      if(logfile) fprintf(logfile,"INFO: more than %d touching neighbor atoms found, growing shear history...",maxtouch);
  }

  int **partner_g = memory->create_2d_int_array(nmax,maxtouch+DELTA_MAXTOUCH,
				      "shear_history:partner_g");
  double ***shearpartner_g =
    memory->create_3d_double_array(nmax,maxtouch+DELTA_MAXTOUCH,3,"shear_history:shearpartner_g");

  for (int i=0;i<nmax;i++)
  {
      for (int j=0;j<maxtouch;j++)
      {
          partner_g[i][j]=partner[i][j];
          for (int k=0;k<3;k++) shearpartner_g[i][j][k]=shearpartner_g[i][j][k];
      }
  }
  maxtouch+=DELTA_MAXTOUCH;
  int **h1; double ***h2;			 ;
  h1 = partner;
  h2 = shearpartner;
  partner = partner_g;
  shearpartner = shearpartner_g;
  memory->destroy_2d_int_array(h1);
  memory->destroy_3d_double_array(h2);
  if(comm->me==0)
  {
    if(screen) fprintf(screen,"done!\n");
    if(logfile) fprintf(logfile,"done!\n");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixShearHistory::copy_arrays(int i, int j)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    shearpartner[j][m][0] = shearpartner[i][m][0];
    shearpartner[j][m][1] = shearpartner[i][m][1];
    shearpartner[j][m][2] = shearpartner[i][m][2];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixShearHistory::set_arrays(int i)
{
  npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = shearpartner[i][n][0];
    buf[m++] = shearpartner[i][n][1];
    buf[m++] = shearpartner[i][n][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixShearHistory::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    shearpartner[nlocal][n][0] = buf[m++];
    shearpartner[nlocal][n][1] = buf[m++];
    shearpartner[nlocal][n][2] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixShearHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 4*npartner[i] + 2;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = shearpartner[i][n][0];
    buf[m++] = shearpartner[i][n][1];
    buf[m++] = shearpartner[i][n][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixShearHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    shearpartner[nlocal][n][0] = extra[nlocal][m++];
    shearpartner[nlocal][n][1] = extra[nlocal][m++];
    shearpartner[nlocal][n][2] = extra[nlocal][m++];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixShearHistory::maxsize_restart()
{
  return 4*maxtouch + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixShearHistory::size_restart(int nlocal)
{
  return 4*npartner[nlocal] + 2;
}
