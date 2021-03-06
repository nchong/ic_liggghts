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

#ifdef PAIR_CLASS

PairStyle(gran/hertz/history/gpu,PairGranHertzHistoryGPU)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_GPU_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_GPU_H

#include "pair_gran_hertz_history.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryGPU : public PairGranHertzHistory {
 public:
  PairGranHertzHistoryGPU(LAMMPS *lmp);
  ~PairGranHertzHistoryGPU();
  void compute(int, int);
  //void settings(int, char **);
  void init_style();
  //double memory_usage();

  enum { ONE_NODE, ONE_GPU, MULTI_GPU };

 private:  
};

}
#endif
#endif
