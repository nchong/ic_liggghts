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

FixStyle(pour/dev/packing,FixPourDevPacking)

#else

#ifndef LMP_FIX_POUR_DEV_PACKING_H
#define LMP_FIX_POUR_DEV_PACKING_H

#include "fix_pour_dev.h"

namespace LAMMPS_NS {

class FixPourDevPacking : public FixPourDev {
  friend class PairGranHertzHistory;
  friend class PairGranHooke;
  friend class PairGranHookeHistory;
  friend class FixPourMultiSphere;
  friend class MechParamGran;

 public:
  FixPourDevPacking(class LAMMPS *, int, char **);
  ~FixPourDevPacking();

  void init();

 private:
  virtual void random_insert_height(double &,double &,double); 
  virtual void calc_insert_velocities(int,double**,double &,double &,double &); 
};

}

#endif
#endif

