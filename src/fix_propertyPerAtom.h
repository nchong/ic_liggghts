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

FixStyle(property/peratom,FixPropertyPerAtom)

#else

#ifndef LMP_FIX_PROPERTYPERATOM_H
#define LMP_FIX_PROPERTYPERATOM_H
#include "fix.h"
#include "input.h"

namespace LAMMPS_NS {

class FixPropertyPerAtom : public Fix {
 public:
  FixPropertyPerAtom(class LAMMPS *, int, char **);
  ~FixPropertyPerAtom();
  int setmask();

  void do_forward_comm();
  void do_reverse_comm();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

  char *variablename;   
  int vectorStyle;      
  int commGhost;        
  int commGhostRev;     
  int nvalues;
  double *defaultvalues; 

 private:

}; //end class

}
#endif
#endif
