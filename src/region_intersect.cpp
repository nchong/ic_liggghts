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

#include "stdlib.h"
#include "string.h"
#include "region_intersect.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

RegIntersect::RegIntersect(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal region command");
  int n = atoi(arg[2]);
  if (n < 2) error->all("Illegal region command");
  options(narg-(n+3),&arg[n+3]);

  // build list of regions to intersect

  list = new int[n];
  nregion = 0;

  int iregion;
  for (int iarg = 0; iarg < n; iarg++) {
    iregion = domain->find_region(arg[iarg+3]);
    if (iregion == -1) error->all("Region intersect region ID does not exist");
    list[nregion++] = iregion;
  }

  // extent of intersection of regions
  // has bounding box if interior and any sub-region has bounding box

  Region **regions = domain->regions;

  bboxflag = 0;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->bboxflag == 1) bboxflag = 1;
  if (!interior) bboxflag = 0;

  if (bboxflag) {
    int first = 1;
    for (int ilist = 0; ilist < nregion; ilist++) {
      if (regions[list[ilist]]->bboxflag == 0) continue;
      if (first) {
	extent_xlo = regions[list[ilist]]->extent_xlo;
	extent_ylo = regions[list[ilist]]->extent_ylo;
	extent_zlo = regions[list[ilist]]->extent_zlo;
	extent_xhi = regions[list[ilist]]->extent_xhi;
	extent_yhi = regions[list[ilist]]->extent_yhi;
	extent_zhi = regions[list[ilist]]->extent_zhi;
	first = 0;
      }

      extent_xlo = MAX(extent_xlo,regions[list[ilist]]->extent_xlo);
      extent_ylo = MAX(extent_ylo,regions[list[ilist]]->extent_ylo);
      extent_zlo = MAX(extent_zlo,regions[list[ilist]]->extent_zlo);
      extent_xhi = MIN(extent_xhi,regions[list[ilist]]->extent_xhi);
      extent_yhi = MIN(extent_yhi,regions[list[ilist]]->extent_yhi);
      extent_zhi = MIN(extent_zhi,regions[list[ilist]]->extent_zhi);
    }
  }

  // possible contacts = sum of possible contacts in all sub-regions

  cmax = 0;
  for (int ilist = 0; ilist < nregion; ilist++)
    cmax += regions[list[ilist]]->cmax;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegIntersect::~RegIntersect()
{
  delete [] list;
  delete [] contact;
}

/* ----------------------------------------------------------------------
   return 1 if region is dynamic, 0 if static
   dynamic if any sub-region is dynamic, else static
------------------------------------------------------------------------- */

int RegIntersect::dynamic_check()
{
  Region **regions = domain->regions;
  for (int ilist = 0; ilist < nregion; ilist++)
    if (regions[list[ilist]]->dynamic_check()) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is match() with all sub-regions
   else inside = 0
------------------------------------------------------------------------- */

int RegIntersect::inside(double x, double y, double z)
{
  int ilist;
  Region **regions = domain->regions;
  for (ilist = 0; ilist < nregion; ilist++)
    if (!regions[list[ilist]]->match(x,y,z)) break;

  if (ilist == nregion) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) compute contacts in each sub-region
   (2) only keep a contact if surface point is match() to all other regions
------------------------------------------------------------------------- */

int RegIntersect::surface_interior(double *x, double cutoff)
{
  int m,ilist,jlist,iregion,jregion,ncontacts;
  double xs,ys,zs;

  Region **regions = domain->regions;
  int n = 0;

  for (ilist = 0; ilist < nregion; ilist++) {
    iregion = list[ilist];
    ncontacts = regions[iregion]->surface(x[0],x[1],x[2],cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = x[0] - regions[iregion]->contact[m].delx;
      ys = x[1] - regions[iregion]->contact[m].dely;
      zs = x[2] - regions[iregion]->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
	if (jlist == ilist) continue;
	jregion = list[jlist];
	if (!regions[jregion]->match(xs,ys,zs)) break;
      }
      if (jlist == nregion) {
	contact[n].r = regions[iregion]->contact[m].r;
	contact[n].delx = regions[iregion]->contact[m].delx;
	contact[n].dely = regions[iregion]->contact[m].dely;
	contact[n].delz = regions[iregion]->contact[m].delz;
	n++;
      }
    }
  }
  
  return n;
}

/* ----------------------------------------------------------------------
   compute contacts with interior of intersection of sub-regions
   (1) flip interior/exterior flag of each sub-region
   (2) compute contacts in each sub-region
   (3) only keep a contact if surface point is not match() to all other regions
   (4) flip interior/exterior flags back to original settings
   this is effectively same algorithm as surface_interior() for RegUnion
------------------------------------------------------------------------- */

int RegIntersect::surface_exterior(double *x, double cutoff)
{
  int m,ilist,jlist,iregion,jregion,ncontacts;
  double xs,ys,zs;

  Region **regions = domain->regions;
  int n = 0;

  for (ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->interior ^= 1;

  for (ilist = 0; ilist < nregion; ilist++) {
    iregion = list[ilist];
    ncontacts = regions[iregion]->surface(x[0],x[1],x[2],cutoff);
    for (m = 0; m < ncontacts; m++) {
      xs = x[0] - regions[iregion]->contact[m].delx;
      ys = x[1] - regions[iregion]->contact[m].dely;
      zs = x[2] - regions[iregion]->contact[m].delz;
      for (jlist = 0; jlist < nregion; jlist++) {
	if (jlist == ilist) continue;
	jregion = list[jlist];
	if (regions[jregion]->match(xs,ys,zs)) break;
      }
      if (jlist == nregion) {
	contact[n].r = regions[iregion]->contact[m].r;
	contact[n].delx = regions[iregion]->contact[m].delx;
	contact[n].dely = regions[iregion]->contact[m].dely;
	contact[n].delz = regions[iregion]->contact[m].delz;
	n++;
      }
    }
  }

  for (ilist = 0; ilist < nregion; ilist++)
    regions[list[ilist]]->interior ^= 1;
  
  return n;
}
