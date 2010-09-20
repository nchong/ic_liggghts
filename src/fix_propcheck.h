/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(propcheck,FixPropCheck)

#else

#ifndef LMP_FIX_PROPCHECK_H
#define LMP_FIX_PROPCHECK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropCheck : public Fix {
 public:
  FixPropCheck(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();

 private:
  int whichpropflag;
  int opflag;
  int actionflag; // not really used..
  double prop_threshold;
  double valtoset;

};

}

#endif
#endif
