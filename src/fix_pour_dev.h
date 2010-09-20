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

#ifdef FIX_CLASS

FixStyle(pour/dev,FixPourDev)

#else

#ifndef LMP_FIX_POUR_DEV_H
#define LMP_FIX_POUR_DEV_H

#include "fix.h"

enum{RAN_STYLE_CONSTANT_FP,RAN_STYLE_UNIFORM_FP,RAN_STYLE_GAUSSIAN_FP}; 

namespace LAMMPS_NS {

class FixPourDev : public Fix {
  friend class PairGranHertzHistory;
  friend class PairGranHooke;
  friend class PairGranHookeHistory;
  friend class FixPourMultiSphere;
  friend class MechParamGran;

 public:
  FixPourDev(class LAMMPS *, int, char **);
  ~FixPourDev();
  int setmask();
  virtual void init();
  void pre_exchange();
  void reset_dt();
  void write_restart(FILE *);
  void restart(char *);
  virtual double max_rad(int);  

 protected:
  class FixParticledistributionDiscrete *fpdd;
  int ninsert;
  int nfinal;
  int nBody;
  int seed;
  int ntype_max; 
  int ntype_min; 
  int check_ol_flag; 
  double masstotal;
  double massflowrate,mass_ins;
  double volfrac;
  int maxattempt;
  int region_style;
  int vel_rand_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;
  double grav;
  int iarg; 

  int me,nprocs;
  int *recvcounts,*displs;
  double PI;
  int nfreq,nfirst,ninserted;
  double nper; 
  double lo_current,hi_current;
  class FixShearHistory *fix_history;
  class RanPark *random;

  //region exempts, 
  int nRegEx;
  class Region **regExList;
  int isInExemptRegion(double *);

  int overlap(int);
  void xyz_random(double, double *);

  double rand_pour(double, double, int); 
  virtual int insert_in_xnear(double **xnear,int nnear,double *coord); 
  virtual bool overlaps_xnear_i(double *coord,double **xnear,int i); 
  virtual int particles_per_insertion(); 
  virtual void random_insert_height(double &,double &,double); 
  virtual void calc_insert_velocities(int,double**,double &,double &,double &); 
  virtual double shift_randompos(); 
  virtual bool give_mol_id(); 
  virtual void init_substyle(){} 
  virtual void set_body_props(bool,int, double,double,double,int){} 
  virtual void finalize_insertion(){} 
};

}

#endif
#endif
