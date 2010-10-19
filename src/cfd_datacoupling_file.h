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

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(file,CfdDatacouplingFile)

#else

#ifndef LMP_CFD_DATACOUPLING_FILE_H
#define LMP_CFD_DATACOUPLING_FILE_H

#include "cfd_datacoupling.h"

namespace LAMMPS_NS {

class CfdDatacouplingFile : public CfdDatacoupling {
 public:
  CfdDatacouplingFile(class LAMMPS *, int, int, char **,class FixCfdCoupling* fc);
  ~CfdDatacouplingFile();

  void pull(char *,char *,void *&);
  void push(char *,char *,void *&);

  private:
   char* filepath;
   bool firstexec;

   char * getFilePath(char *name,bool);
   void op_complete(char *name);
   void writeVectorData(char *name,  double ** field);
   void writeScalarData(char *name,  double * field);
   void writeGlobalVectorData(char *name,  double * field,int);
   void writeGlobalArrayData(char *name,  double ** field,int,int);

   void readVectorData(char *name,  double ** field);
   void readScalarData(char *name,  double * field);
   void readGlobalVectorData(char* name, double *field, int &len);
   void readGlobalArrayData(char *name, double ** field, int &len1, int &len2);

};

}

#endif
#endif
