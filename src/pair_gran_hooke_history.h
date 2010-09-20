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

#ifdef PAIR_CLASS

PairStyle(gran/hooke/history,PairGranHookeHistory)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_H

#include "pair.h"
#include "fix_propertyGlobal.h"

namespace LAMMPS_NS {

class PairGranHookeHistory : public Pair {
 public:
  friend class FixWallGranHookeHistory;
  friend class FixCheckTimestepGran;
  PairGranHookeHistory(class LAMMPS *);
  ~PairGranHookeHistory();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual void init_substyle(); 
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void reset_dt();

 protected:

   class MechParamGran *mpg;

   class FixRigid* fr;

  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &);
  virtual void addCohesionForce(int &, int &,double &,double &);

  int cohesionflag; 

  virtual class MechParamGran* getMatProp();
  
  int dampflag;
  double dt;
  int freeze_group_bit;
  int history;
  class FixShearHistory *fix_history;

  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  void allocate();
};

}

#endif
#endif
