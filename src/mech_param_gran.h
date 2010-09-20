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

#ifndef LMP_MECHPARAMGRAN_H
#define LMP_MECHPARAMGRAN_H

#include "lammps.h"
#include "pointers.h"

namespace LAMMPS_NS {

class MechParamGran: protected Pointers
{
 friend class PairGranHookeHistory;
 friend class FixWallGranHookeHistory;
 friend class FixCheckTimestepGran;

 public:

  MechParamGran(LAMMPS *lmp);
  ~MechParamGran();

  int min_type,max_type;

  class FixPropertyGlobal* Y1; //Youngs Modulus
  class FixPropertyGlobal* v1; //Poisson's ratio
  class FixPropertyGlobal* cohEnergyDens1; //Cohesion energy density

  class FixPropertyGlobal* coeffRest1; //coefficient of restitution
  class FixPropertyGlobal* coeffFrict1; //coefficient of (static) friction

  class FixPropertyGlobal* charVel1; //characteristic velocity needed for Linear Spring Model

  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,charVel;

  bool arrays_active;

  void create_arrays(int);
  void destroy_arrays();
  void getMaterialParams(int,int);
}; //end class

}

#endif
