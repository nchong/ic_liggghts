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

#ifndef LMP_CFD_REGIONMODEL_H
#define LMP_CFD_REGIONMODEL_H

#include "pointers.h"

namespace LAMMPS_NS {

class CfdRegionmodel : protected Pointers {
 public:
  CfdRegionmodel(class LAMMPS *lmp, int jarg,int narg, char **arg,class FixCfdCoupling* fc) : Pointers(lmp)
  {
      this->fc = fc;
  }
  ~CfdRegionmodel() {}

  int get_iarg() {return iarg;}
  bool liggghts_is_active;

  virtual void init() {};
  virtual void rm_update() {};

 protected:
  int iarg;
  class FixCfdCoupling *fc;

  void add_push_property(char *name,char *type)
  {
     fc->add_push_property(name,type);
  }

  int add_pull_property(char *name,char *type)
  {
     fc->add_pull_property(name,type);
  }
};

}

#endif

