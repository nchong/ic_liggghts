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

FixStyle(particletemplate/sphere,FixTemplateSphere)

#else

#ifndef LMP_FIX_TEMPLATE_SPHERE_H
#define LMP_FIX_TEMPLATE_SPHERE_H

#include "fix.h"

enum{RAN_STYLE_CONSTANT_FTS,RAN_STYLE_UNIFORM_FTS,RAN_STYLE_GAUSSIAN_FTS};

namespace LAMMPS_NS {

class FixTemplateSphere : public Fix {
 public:
  FixTemplateSphere(class LAMMPS *, int, char **);
  ~FixTemplateSphere();
  void write_restart(FILE *);
  void restart(char *);

  virtual int setmask();
  virtual double volexpect();           
  virtual double massexpect();          
  double max_rad();

  virtual void randomize();              
  class ParticleToInsert* pti;

 protected:
  class RanPark *random;
  int seed;
  double PI;
  int iarg;

  //properties of particle template
  int nspheres;       
  double **x_sphere;   
  double **x_sphere_b; 
  double *r_sphere;   

  int atom_type;
  double density1;          // = min val for uniform, µ for gauss
  double density2;          // = max val for uniform, sigma for gauss
  int density_randstyle;
  double volume;
  double mass;
  int radius_randstyle;
  double radius2;               //= ratio max/min for uniform, sigma for gauss

  virtual void randomize_r();
  virtual void randomize_dens();

};

}

#endif
#endif

