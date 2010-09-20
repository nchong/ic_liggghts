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

#include "string.h"
#include "dump_stl.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_meshGran.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define nFixMeshGran 20 //a maximum of nFixMeshGran fixes of type 'mesh/gran' is allowed

/* ---------------------------------------------------------------------- */

DumpSTL::DumpSTL(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal dump stl command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  //if (multifile!=1) error->warning("You should use a filename like 'dump*.stl' for the 'dump stl' command if you want one file per time-step");
  if (multiproc!=0) error->warning("Your 'dump stl' command is writing one file per processor, where all the files contain the same data");

  STLList=new STLtri*[nFixMeshGran];
  for (int i=0;i<nFixMeshGran;i++) STLList[i]=(STLtri*)NULL;

  format_default = NULL;

  size_one=12;  

  int iarg=5;
  meshid_len=narg-iarg;
  meshid=new char*[meshid_len];
  for (int k=0;k<meshid_len;k++) meshid[k]=new char[40];

  int k=0;
  while(iarg<narg)
      strcpy(meshid[k++],arg[iarg++]);

}

/* ---------------------------------------------------------------------- */

DumpSTL::~DumpSTL()
{
    for (int k=0;k<meshid_len;k++) delete []meshid[k];
}

/* ---------------------------------------------------------------------- */

void DumpSTL::init()
{
  if (domain->triclinic == 1) error->all("Can not dump STL files for triclinic box");
  if (binary) error->all("Can not dump STL files in binary mode");

  STLList_len=0;
  for (int ifix = 0; ifix < modify->nfix; ifix++)
  {
    if (strncmp(modify->fix[ifix]->style,"mesh/gran",8) == 0)  {

        bool found = false;
        for (int k=0;k<meshid_len;k++)
            if(strcmp(modify->fix[ifix]->id,meshid[k])==0) found=true;
        if(!found && meshid_len > 0) continue;

        STLList[STLList_len]=static_cast<FixMeshGran*>(modify->fix[ifix])->STLdata;
        STLList_len++;
    }
    if (STLList_len==19) error->all("Dump stl can not process all fixes of type mesh/gran, boost nFixMeshGran in dump_stl.cpp and recompile");
  }

  if(STLList_len<meshid_len) error->warning("Dump stl could not locate all fixes of type 'mesh/gran' that you wish to dump");
  if (STLList_len==0) error->warning("Dump stl cannot find any fix of type 'mesh/gran' to dump");

  // default format depends on image flags

  delete [] format;
  format = new char[150];
  strcpy(format,"facet normal %g %g %g\n");
  strcat(format,"outer loop\n");
  strcat(format,"vertex %g %g %g\n");
  strcat(format,"vertex %g %g %g\n");
  strcat(format,"vertex %g %g %g\n");
  strcat(format,"endloop\n");
  strcat(format,"endfacet\n");

  size_one=12; //12 double values per tri

  // setup function ptrs

  header_choice = &DumpSTL::header_item;
  pack_choice = &DumpSTL::pack_item;
  write_choice = &DumpSTL::write_item;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpSTL::modify_param(int narg, char **arg)
{
  error->warning("dump_modify keyword is not supported by 'dump stl' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpSTL::write_header(int ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpSTL::count()
{
  if (comm->me!=0) return 0;
  int nTri=0;
  for(int i=0;i<STLList_len;i++) nTri+=STLList[i]->nTri;
  return nTri;
}

/* ---------------------------------------------------------------------- */

int DumpSTL::pack()
{
  return (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpSTL::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpSTL::header_item(int ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"solid LIGGGHTS_STL_EXPORT\n");
}

void DumpSTL::footer_item()
{
  if (comm->me!=0) return;
  fprintf(fp,"endsolid LIGGGHTS_STL_EXPORT\n");
}

/* ---------------------------------------------------------------------- */

int DumpSTL::pack_item()
{
  if (comm->me!=0) return 0;
    int m = 0;
    STLtri *STLdat;
    for(int i=0;i<STLList_len;i++)
    {
        STLdat=STLList[i];
        int nTri=STLdat->nTri;
        for(int j=0;j<nTri;j++)
        {
            //process triangle j of STLdat
            for(int k=0;k<3;k++) buf[m++]=STLdat->facenormal[j][k];
            for(int k=0;k<3;k++) { //loop all nodes
                for(int l=0;l<3;l++){ //loop x,y,z of a node
                     buf[m++]=STLdat->node[j][k][l];
                }
            }
        }
    }
    return m;
}

/* ---------------------------------------------------------------------- */

void DumpSTL::write_item(int n, double *mybuf) 
{
  if (comm->me!=0) return;
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
	    mybuf[m],mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],
	    mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9],mybuf[m+10],mybuf[m+11]
    );
    m += size_one;
  }
  //write footer
  footer_item();

}
