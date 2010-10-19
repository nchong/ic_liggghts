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

#ifndef LMP_CFD_DATACOUPLING_H
#define LMP_CFD_DATACOUPLING_H

#include "pointers.h"

namespace LAMMPS_NS {

class CfdDatacoupling : protected Pointers {
 public:

  CfdDatacoupling(class LAMMPS *lmp, int jarg,int narg, char **arg,class FixCfdCoupling* fc) : Pointers(lmp)
  {
      this->fc = fc;
  }
  ~CfdDatacoupling() {}

  int get_iarg() {return iarg;}
  bool liggghts_is_active;

  virtual void pull(char *,char *,void *&) {};
  virtual void push(char *,char *,void *&) {};

 protected:
  int iarg;
  class FixCfdCoupling *fc;
};

}

#endif
