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

FixStyle(neighlist/tri,FixTriNeighlist)

#else

#ifndef LMP_FIX_TRI_NEIGHLIST_H
#define LMP_FIX_TRI_NEIGHLIST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTriNeighlist : public Fix {

  friend class FixWallGranHookeHistory;

 public:
  FixTriNeighlist(class LAMMPS *, int, char **);
  ~FixTriNeighlist();
  int setmask();
  void init();
  void pre_neighbor(); 
  void pre_force(int); 

  double memory_usage();
  void grow_arrays(int);
  void grow_arrays_maxtritouch(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 protected:
  
  int *nTriList;

  int ***tri_neighlist;
  int maxwalllist;
  int *delflag;

 private:
  int addTriToNeighList(int,int,int);
  int check_dangerous(double,int);
  void flag_old_list();
  void clear_old_entries();

  void check_dangerous_build(int);
  int do_warn,do_warn_dangerous;

  bool check_box_overlap(double*,double*,double*,double*);

  int decide_rebuild();
  int check_distance();
  void unset_nontouching();

  char *caller_id;
  class FixWallGranHookeHistory* caller;
  int nFixMeshGran;
  class FixMeshGran** FixMeshGranList;
  int buildNeighList;  

  //neighbor list params
  double *bsubboxlo, *bsubboxhi;
};

}

#endif
#endif
