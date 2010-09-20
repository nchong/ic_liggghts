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

FixStyle(mesh/gran/stressanalysis,FixMeshGranAnalyze)

#else

#ifndef LMP_FIX_MESHGRAN_ANALYZE_H
#define LMP_FIX_MESHGRAN_ANALYZE_H

#include "fix.h"
#include "math.h"
#include "fix_meshGran.h"

namespace LAMMPS_NS {

class FixMeshGranAnalyze : public FixMeshGran {

 public:
  FixMeshGranAnalyze(class LAMMPS *, int, char **);
  ~FixMeshGranAnalyze();
  virtual int setmask();
  void pre_force(int);
  void add_particle_contribution(double*,double*,int);
  virtual void final_integrate();
  virtual void init(){};
  virtual int write_restart_sub(FILE * fp,int n){return n;}
  virtual void restart_sub(char *){}
  double compute_vector(int);

 protected:
  void calc_total_force();
  virtual int n_children(){return 0;}
  virtual void children_write(FILE* fp) {}
  virtual void children_restart(double *){}

 private:
  double *tmp,*tmp2;

}; //end class

}

#endif
#endif

