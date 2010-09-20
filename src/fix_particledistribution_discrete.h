
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

FixStyle(particledistribution/discrete,FixParticledistributionDiscrete)

#else

#ifndef LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H
#define LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H

#include "fix.h"

enum{RAN_STYLE_CONSTANT_FPD,RAN_STYLE_UNIFORM_FPD,RAN_STYLE_GAUSSIAN_FPD};

namespace LAMMPS_NS {

class FixParticledistributionDiscrete : public Fix {
 public:
  friend class FixPourDev;
  FixParticledistributionDiscrete(class LAMMPS *, int, char **);
  ~FixParticledistributionDiscrete();
  void write_restart(FILE *);
  void restart(char *);

  int setmask();
  double vol_expect();
  double mass_expect();
  int max_type();
  int min_type();
  double max_rad(int);
  int max_nspheres();

  int ninserted;
  int ninsert;
  int random_init(int);         
  void randomize();             
  class ParticleToInsert* pti;

 protected:
  class RanPark *random;
  int seed;

  int iarg;

  //particle templates
  int ntemplates;
  double *distweight;
  double *cumweight;
  int *parttogen;
  int *distorder;
  class FixTemplateSphere **templates;

  //mass and volume expectancy of this discrete distribution
  double volexpect;
  double massexpect;

  //min/maximum particle type to be inserted
  int maxtype;
  int mintype;

  //maximum number of spheres a template has
  int maxnspheres;
};

}

#endif
#endif

