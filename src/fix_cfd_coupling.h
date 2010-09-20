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

#else

#ifndef LMP_FIX_CFD_COUPLING_H
#define LMP_FIX_CFD_COUPLING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCfdCoupling : public Fix {
 public:
  FixCfdCoupling(class LAMMPS *, int, char **);
  ~FixCfdCoupling();

  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void min_setup(int);
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void pre_force(int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  void end_of_step();
  void post_force_respa(int, int, int);
  void min_post_force(int);

  double compute_array(int,int);

  double memory_usage();
  void grow_arrays(int);
  void grow_arrays_allred(int);
  void zero_arrays();
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void set_arrays(int);

  void pull(char *name,char *type,void *&ptr);
  void push(char *name,char *type,void *&ptr);
  void *find_pull_property(char *name,char *type);

 protected:
  void add_push_property(char *name,char *type);
  int add_pull_property(char *name,char *type);

  virtual void special_settings() =0;

  int nvalues;

 private:

  int couple_nevery;

  int nvalues_max;
  
  int *nvals_each;
  char **valnames;
  char **valtypes;
  
  int npush;
  char **pushnames;
  char **pushtypes;

  void grow_();

  class CfdDatacoupling *dc;

  int nlevels_respa;
};

}

#endif
#endif
