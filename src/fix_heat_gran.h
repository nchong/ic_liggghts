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

FixStyle(heat/gran,FixHeatGran)

#else

#ifndef LMP_FIX_HEATGRAN_H
#define LMP_FIX_HEATGRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeatGran : public Fix {
 public:
  FixHeatGran(class LAMMPS *, int, char **);
  ~FixHeatGran();
  int setmask();
  void init();
  void updatePtrs();
  void post_force(int);
  void final_integrate();
  void initial_integrate_respa(int,int,int);
  void initial_integrate(int);
  double compute_scalar();
  void reset_dt();

 private:
  class FixPropertyPerAtom* fppat;
  class FixPropertyPerAtom* fppahf;
  class FixPropertyPerAtom* fppahs;
  class FixPropertyGlobal* fpgco;
  class FixPropertyGlobal* fpgca;
  double T0;              
  double *Temp;           
  double *heatFlux;       
  double *heatSource;     

};

}

#endif
#endif

