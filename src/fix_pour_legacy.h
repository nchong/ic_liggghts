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

#ifdef FIX_CLASS

FixStyle(pour/legacy,FixPourLegacy)

#else

#ifndef LMP_FIX_POUR_LEGACY_H
#define LMP_FIX_POUR_LEGACY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPourLegacy : public Fix {
  friend class PairGranHertzHistory;
  friend class PairGranHooke;
  friend class PairGranHookeHistory;

 public:
  FixPourLegacy(class LAMMPS *, int, char **);
  ~FixPourLegacy();
  int setmask();
  void init();
  void pre_exchange();
  void reset_dt();
  virtual double max_rad(int);  

 private:
  int ninsert,ntype,seed;
  double radius_lo,radius_hi;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;

  int me,nprocs;
  int *recvcounts,*displs;
  double PI;
  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  class FixShearHistory *fix_history;
  class RanPark *random;

  int overlap(int);
  void xyz_random(double, double *);
};

}

#endif
#endif
