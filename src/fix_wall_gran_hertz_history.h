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

FixStyle(wall/gran/hertz/history,FixWallGranHertzHistory)

#else

#ifndef LMP_FIX_WALL_GRAN_HERTZ_HISTORY_H
#define LMP_FIX_WALL_GRAN_HERTZ_HISTORY_H

#include "fix.h"
#include "fix_wall_gran_hooke_history.h"

namespace LAMMPS_NS {

class FixWallGranHertzHistory : public FixWallGranHookeHistory {
 public:
  FixWallGranHertzHistory(class LAMMPS *, int, char **);
  int setmask();
  void init();

 protected:

  virtual void deriveContactModelParams(int &, double , double &, double &, double &, double &, double &);
  virtual void initSubstyle();
};

}

#endif
#endif
