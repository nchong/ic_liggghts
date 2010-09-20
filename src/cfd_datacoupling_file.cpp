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

#include "sys/stat.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling.h"
#include "cfd_datacoupling_file.h"
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

CfdDatacouplingFile::CfdDatacouplingFile(LAMMPS *lmp, int narg, char **arg,FixCfdCoupling *fc)  :
  CfdDatacoupling(lmp, narg, arg,fc)
{
    if(comm->nprocs>1)  error->all("Fix couple/cfd with file coupling is for serial computation only");
    if(narg==0) error->all("Fix couple/cfd: not enough arguments");

    liggghts_is_active = true;
    firstexec = true;

    this->fc = fc;

    filepath = new char[strlen(arg[0])+2];
    strcpy(filepath,arg[0]);
    if(filepath[strlen(arg[0])]!='/') strcat(filepath,"/");

}

CfdDatacouplingFile::~CfdDatacouplingFile()
{
    delete []filepath;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::pull(char *name,char *type,void *&from)
{
    void * to = fc->find_pull_property(name,type);

    if (!to)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->all("This is fatal");
    }

    if(strcmp(type,"scalar") == 0)
    {
        readScalarData(name,(double*)to);
    }

    if(strcmp(type,"vector") == 0)
    {
        readVectorData(name,(double**)to);
    }
    firstexec = false;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingFile::push(char *name,char *type,void *&to)
{
    void *from = atom->extract(name);

    if (!from)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s requested by calling program.\n",name);
        lmp->error->all("This is fatal");
    }

    if(strcmp(type,"scalar") == 0)
    {
        writeScalarData(name,(double*)from);
    }

    if(strcmp(type,"vector") == 0)
    {
        writeVectorData(name,(double**)from);
    }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

char * CfdDatacouplingFile::getFilePath(char *name,bool flag)
{
    char *file = new char[strlen(filepath)+strlen(name)+3];
    strcpy(file,filepath);
    strcat(file,name);
    if(flag) strcat(file,"0");
    else strcat(file,"1");
    struct stat st;
    return file;
}

void CfdDatacouplingFile::op_complete(char *name)
{
    char *oldfile = getFilePath(name,true);
    char *newfile = getFilePath(name,false);
    
    rename(oldfile,newfile);
    delete []oldfile;
    delete []newfile;
}

void CfdDatacouplingFile::readVectorData(char *name, double ** field)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(0.03);

    // set file pointer
    ifstream inputPtr(file);

    // write data to variable
    int numberOfParticles;
    inputPtr >> numberOfParticles;
    
    if(atom->nlocal!=numberOfParticles) error->all("Fix couple/cfd/file: Data corruption: # particles in file does not match # particles in LIGGGHTS");

    for(int index = 0;index < numberOfParticles; ++index)
    {
        for(int i=0;i<3;i++) inputPtr >> field[index][i];
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

void CfdDatacouplingFile::readScalarData(char* name, double *field)
{
    // get output path
    char *file = getFilePath(name,true);

    fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
    struct stat st;
    while (stat(file,&st)) sleep(0.03);

    // set file pointer
    ifstream inputPtr(file);

    // write data to variable
    int numberOfParticles;
    inputPtr >> numberOfParticles;

    if(atom->nlocal!=numberOfParticles) error->all("Fix couple/cfd/file: Data corruption: # particles in file does not match # particles in LIGGGHTS");

    // write data to variable
    for(int index = 0;index < numberOfParticles; ++index)
    {
        inputPtr >> field[index];
    }

    // clean up inputStream
    delete []file;
    op_complete(name);
}

void CfdDatacouplingFile::writeVectorData(char *name,  double ** field)
{
    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(0.03);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    int numberOfParticles = atom->nlocal;
    outputPtr << numberOfParticles << endl;
    for(int index = 0;index < numberOfParticles; ++index)
    {
        for(int i=0;i<3;i++) outputPtr << field[index][i] << " ";
        outputPtr << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}

void CfdDatacouplingFile::writeScalarData(char* name, double * field)
{
    // get output path
    char *file = getFilePath(name,true);

    if(!firstexec)
    {
      fprintf(screen,"Fix couple/cfd/file: waiting for file: %s\n",file);
       struct stat st;
       while (stat(file,&st)) sleep(0.03);
    }

    // set file pointer
    ofstream outputPtr(file);

    // write data to file
    int numberOfParticles = atom->nlocal;
    outputPtr << numberOfParticles << endl;
    for(int index = 0;index < numberOfParticles; ++index)
    {
        outputPtr << field[index] << endl;
    }

    // clean up outputStream and rename file
    delete []file;
    op_complete(name);
}

