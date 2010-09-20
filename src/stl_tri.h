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

#ifndef LMP_STLTRI_H
#define LMP_STLTRI_H

#include "lammps.h"
#include "pointers.h"

namespace LAMMPS_NS {

class STLtri: protected Pointers {

 public:
  STLtri(LAMMPS *lmp);
  ~STLtri();

  void grow_arrays();
  void pack_restart();
  void unpack_restart();
  void initMove(double);
  void initConveyor();
  void pointToVec();
  void vecToPoint();
  void getTriBoundBox(int, double *, double *,double);
  void before_rebuild();
  bool are_coplanar_neighs(int i1,int i2);

  int nTri,nTriMax;

  int movingMesh;  
  int conveyor;    
  double conv_vel[3];
  double skinSafetyFactor; 

  int* EDGE_INACTIVE;
  int* CORNER_INACTIVE;

   #define VECPERTRI 11

   //these are the triangles read from the STL file
   //index1:facet(triangle) index2:x y z

   double ***node;        
   double ***v_node;      
   double ***node_lastRe; 
   double **facenormal;
   double **f_tri;        
   double **fn_fshear;    
   double *Area;          

   double **cK;
   double ***ogK;
   double **ogKlen;
   
   double ***oKO; 
   double *rK;

   int **neighfaces; 
   int  *contactInactive;

  //the data is stored in the arrays above, the x array below have pointers to the data
  
  double **x;         
  double **xoriginal; 
  
  double **v;
  double **f;
  double *rmass;  
  int xvf_len;    
  int xvf_lenMax; 

 private:
  bool alreadyInitialized;

}; //end class

}

#endif

