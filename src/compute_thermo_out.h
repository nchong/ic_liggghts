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

#ifdef COMPUTE_CLASS

ComputeStyle(thermo/out,ComputeThermoOut)

#else

#ifndef LMP_COMPUTE_THERMO_OUT_H
#define LMP_COMPUTE_THERMO_OUT_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeThermoOut : public Compute {
 public:
  ComputeThermoOut(class LAMMPS *, int, char **);
  ~ComputeThermoOut();
  virtual void init();
  void callback(int ifield,bool isint,int intval,double doubleval,char *keyword);

 private:

  int me;
  FILE *fp,*fp_h;
  char *filename, *filename_h;
  int thermo_vals,called_first,called_last;
  bool is_first_call;
};

}

#endif
#endif
