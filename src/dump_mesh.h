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

#ifdef DUMP_CLASS

DumpStyle(mesh/gran/VTK,DumpMesh)

#else

#ifndef LMP_DUMP_MESH_H
#define LMP_DUMP_MESH_H

#include "dump.h"
#include "stl_tri.h"

enum{DUMP_STRESS,DUMP_ID};

namespace LAMMPS_NS {

class DumpMesh : public Dump {
 public:
  DumpMesh(LAMMPS *, int, char**);
  ~DumpMesh();
  void init();

 private:
  STLtri** STLList;
  int STLList_len;
  char **meshid;
  int meshid_len;
  int dump_what;

  char *columns;             // column labels
  int lasttimestep;

  int modify_param(int, char **);
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  typedef void (DumpMesh::*FnPtrHeader)(int);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_item(int);
  void footer_item();

  typedef int (DumpMesh::*FnPtrPack)();
  FnPtrPack pack_choice;               // ptr to pack functions
  int pack_item();

  typedef void (DumpMesh::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_item(int, double *);

};

}

#endif
#endif
