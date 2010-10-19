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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/hooke/history,FixWallGranHookeHistory)

#else

#ifndef LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H
#define LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H

#include "fix.h"

enum{MESHGRAN_FWGHH,XPLANE_FWGHH,YPLANE_FWGHH,ZPLANE_FWGHH,ZCYLINDER_FWGHH};

#define F_SHRINKAGE -0.000000001  

namespace LAMMPS_NS {

class FixWallGranHookeHistory : public Fix {
  friend class FixTriNeighlist; 
  friend class MechParamGran; 
  friend class FixCheckTimestepGran; 
  friend class FixConstrainMeshGran6DOF; 

 public:
  FixWallGranHookeHistory(class LAMMPS *, int, char **);
  ~FixWallGranHookeHistory();
  virtual int setmask();
  virtual void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);

  double memory_usage();
  void grow_arrays(int);
  virtual void grow_arrays_maxtritouch(int); 
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  void write_restart(FILE *);
  void restart(char *);      
  int size_restart(int);
  int maxsize_restart();
  void reset_dt();

 protected:
  int atom_type_wall;
  double Temp_wall;
  int wallstyle,pairstyle,wiggle,wshear,axis;
  int dampflag,cohesionflag;

  int nFixMeshGran;
  class FixMeshGran** FixMeshGranList;

  double lo,hi,cylradius;
  double amplitude,period,omega,vshear;
  double dt;
  int nlevels_respa;
  int time_origin;

  //list of contact partners
  int maxpartners; 
  int *npartners;
  int ***partner;

  //shear history
  int shearhistory;
  double ***shear;

  class FixRigid *fr;

  class FixTriNeighlist* fix_tri_neighlist;

  class MechParamGran *mpg;

  class FixPropertyPerAtom *fppa_T;
  class FixPropertyPerAtom *fppa_hf;
  double *Temp_p;
  double *heatflux;
  const double *th_cond;
  double const* const* deltan_ratio;

  void init_heattransfer();
  void addHeatFlux(int, double,double, int);
  virtual void compute_force(int, double, double, double, double, double *, double *, double *, double *, double *, double, double, double *,
                                  double, double, double, double, double);
  virtual void addCohesionForce(int &, double &, double &);
  virtual void deriveContactModelParams(int &, double , double &, double &, double &, double &, double &);
  virtual void resetShearHistory(int,int);
  virtual void initSubstyle();
  void registerTriNeighlist(int);

  void reset_wall_forces();
  virtual int add_to_contact_list(int, int, int);
  void shear_transition(int,int,int);
  virtual void remove_from_contact_list(int, int, int);
};

}

#endif
#endif
