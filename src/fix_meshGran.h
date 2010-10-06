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

FixStyle(mesh/gran,FixMeshGran)

#else

#ifndef LMP_FIX_MESHGRAN_H
#define LMP_FIX_MESHGRAN_H

#include "fix.h"
#include "input.h"
#include "math.h"
#include "myvector.h"
#include "stl_tri.h"

namespace LAMMPS_NS {

class FixMeshGran : public Fix {
  friend class FixMoveTri;
  friend class FixTriNeighlist;
  friend class Input;
  friend class DumpSTL;
  friend class DumpMesh;
  friend class FixWallGranHookeHistory;
  friend class FixWallGranHooke;
  friend class FixWallGranHertzHistory;

 public:
  FixMeshGran(class LAMMPS *, int, char **);
  ~FixMeshGran();
  virtual int setmask();
  void write_restart(FILE *);
  void restart(char *);
  virtual int write_restart_sub(FILE * fp,int n){return n;}
  virtual void restart_sub(char *){}

  virtual void final_integrate(){}
  virtual void init(){};
  virtual void post_integrate(){}

  int atom_type_wall;//atom type that is assigned to the wall (to derive the mechanical properties) and thus to all pseudo wall particles

  //these are the triangles read from the STL file
  //index1:facet(triangle) index2:x y z
  class STLtri *STLdata;
  double ***node;
  int nTri;

  double *force_total;
  double *torque_total;
  double *p_ref;

  virtual void add_particle_contribution(double*,double*,int){}

 protected:
  bool analyseStress;
  int iarg;
  double scale_fact,*off_fact, *rot_angle; 

  virtual int n_children(){return 0;}
  virtual void children_write(FILE* fp) {}
  virtual void children_restart(double *){}

 private:
  
  int* EDGE_INACTIVE;
  int* CORNER_INACTIVE;

  char* stl_filename;
  FILE* stl_file;
  class Input *mystl_input;
  void calcTriCharacteristics(int nTri,double ***node,double **cK,double ***ogK,double **ogKlen,double ***oKO,double *rK,double *Area,double **facenormal,int **neighfaces,int *contactInactive);
  int get_max_index_sharedcorner(int iTri,int &nPrev,int *prevFaces,double *node2check,double ***node,double *rK,int **neighfaces);
}; //end class

}

#endif
#endif
