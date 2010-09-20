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

#ifndef LMP_MYVECTOR_H
#define LMP_MYVECTOR_H

#include<cmath>

namespace LAMMPS_NS {
//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

inline double vectorMag3D(double *v)
{
  return (  std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])  );
}

inline double vectorMag3DSquared(double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]  );
}

inline double vectorDot3D(double *v1,double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

inline void vectorCopy3D(double *from,double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorScalarMult3D(double *v,double s)
{
  v[0]=s*v[0];
  v[1]=s*v[1];
  v[2]=s*v[2];
}

inline void vectorScalarMult3D(double *v,double s,double *v2)
{
  v2[0]=s*v[0];
  v2[1]=s*v[1];
  v2[2]=s*v[2];
}

inline void vectorScalarDiv3D(double *v,double s)
{
  v[0]=1./s*v[0];
  v[1]=1./s*v[1];
  v[2]=1./s*v[2];
}

inline void vectorScalarDiv3D(double *v,double s,double *v2)
{
  v2[0]=1./s*v[0];
  v2[1]=1./s*v[1];
  v2[2]=1./s*v[2];
}

inline void vectorAdd3D(const double *v1,const double *v2, double *v3)
{
  v3[0]=v1[0]+v2[0];
  v3[1]=v1[1]+v2[1];
  v3[2]=v1[2]+v2[2];
}

inline void vectorSubtract3D(const double *v1,const double *v2, double *v3)
{
  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];
  v3[2]=v1[2]-v2[2];
}

inline void vectorCross3D(const double *v1,const double *v2, double *v3)
{
  v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

inline void normalize_bary(double *v)
{
  double mag = v[0]+v[1]+v[2];
  v[0]/=mag;
  v[1]/=mag;
  v[2]/=mag;
}

}

#endif
