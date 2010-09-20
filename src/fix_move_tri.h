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

FixStyle(move/mesh/gran,FixMoveTri)

#else

#ifndef LMP_FIX_MOVE_TRI_H
#define LMP_FIX_MOVE_TRI_H

#include "fix.h"
#include "stl_tri.h"

namespace LAMMPS_NS {

class FixMoveTri : public Fix {
 public:
  FixMoveTri(class LAMMPS *, int, char **);
  ~FixMoveTri();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int);
  void write_restart(FILE *);
  void restart(char *);

 private:
  char *xvarstr,*yvarstr,*zvarstr,*vxvarstr,*vyvarstr,*vzvarstr;
  int mstyle;
  int vxflag,vyflag,vzflag,axflag,ayflag,azflag;
  double vx,vy,vz,ax,ay,az;
  double period,omega_rotate,ampl; 
  double point[3],axis[3],runit[3];
  double dt,dtv,dtf;
  int xvar,yvar,zvar,vxvar,vyvar,vzvar;
  int xvarstyle,yvarstyle,zvarstyle,vxvarstyle,vyvarstyle,vzvarstyle;
  int omega_flag,nlevels_respa;
  int time_origin;
  double skinSafetyFactor;            
  class STLtri* STL_tri;              
  double **x,**v,**xoriginal,**f;     
  double *rmass,*mass;                
  int nlocal, *type;                  
  int displaceflag,velocityflag;
  int maxatom;
  double **displace,**velocity;
};

}

#endif
#endif
