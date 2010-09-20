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

FixStyle(pour,FixPour)

#else

#ifndef LMP_FIX_POUR_H
#define LMP_FIX_POUR_H

#include "fix.h"

enum{RAN_STYLE_UNIFORM,RAN_STYLE_GAUSSIAN,RAN_STYLE_LOGNORMAL,RAN_STYLE_DISCRETE}; 

namespace LAMMPS_NS {

class FixPour : public Fix {
  friend class PairGranHertzHistory;
  friend class PairGranHooke;
  friend class PairGranHookeHistory;
  friend class FixPourMultiSphere;
  friend class MechParamGran;

 public:
  FixPour(class LAMMPS *, int, char **);
  ~FixPour();
  int setmask();
  void init();
  virtual void init_substyle() {}; 
  void pre_exchange();
  void reset_dt();
  virtual double max_rad(int); 

 protected:
  virtual void finalize_insertion(){}

 private:
  int ninsert,ntype,seed;
  int radius_ran_style,density_ran_style,vel_ran_style;
  double radius_lo,radius_hi;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;
  double grav; //gravity in -z direction 

  int me,nprocs;
  int *recvcounts,*displs;
  double PI;
  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  class FixShearHistory *fix_history;
  class RanPark *random;

  int overlap(int);
  void xyz_random(double, double *,double);
  virtual void calc_nper();  
  virtual void calc_nfreq(); 
  double rand_pour(double, double, int); 
  double expectancy(double,double,int); 
  double volume_expectancy(double,double,int,int); 
  virtual double density_scaling(); 
  virtual int insert_in_xnear(double **xnear,int nnear,double *coord,double radtmp); 
  virtual bool overlaps_xnear_i(double *coord,double radtmp,double **xnear,int i); 
  virtual int particles_per_insertion(); 
  virtual void calc_insert_velocities(int,double,double**,double &,double &,double &); 
  virtual void set_body_props(int,int,double*,double,double,double,double,double,int){};

  //particles shall be completely within insertion volume
  //possible random positions are shifted depending on the particle radius,
  virtual double shift_randompos(double);  
  void check_gravity(); 
};

}

#endif
#endif
