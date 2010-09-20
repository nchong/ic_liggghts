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

ComputeStyle(totalmass,ComputeTotalMass)

#else

#ifndef LMP_COMPUTE_TOTALMASS_H
#define LMP_COMPUTE_TOTALMASS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTotalMass : public Compute {
 public:
  ComputeTotalMass(class LAMMPS *, int, char **);
  ~ComputeTotalMass();
  void init();
  virtual double compute_scalar();
  virtual void compute_local();

 private:
 int nevery,me;
 FILE *fp;
 char* filename;
 double time;
  //double masstotal; stored in scalar
};

}

#endif
#endif
