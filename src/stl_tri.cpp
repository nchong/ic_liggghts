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

#include "stl_tri.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"
#include "math.h"
#include "domain.h"
#include "myvector.h"

#define DELTATRI 5000 

using namespace LAMMPS_NS;

STLtri::STLtri(LAMMPS *lmp): Pointers(lmp)
{
    node=NULL;
    node_lastRe=NULL;
    v_node=NULL;
    f_tri=NULL;
    fn_fshear=NULL;
    facenormal=NULL;
    cK=NULL;
    ogK=NULL;
    ogKlen=NULL;
    oKO=NULL;
    rK=NULL;
    x=NULL;
    v=NULL;
    xoriginal=NULL;
    f=NULL;
    rmass=NULL;
    Area=NULL;
    neighfaces=NULL;
    contactInactive=NULL;

    nTri=0;
    nTriMax=0;
    xvf_len=0;
    xvf_lenMax=0;
    grow_arrays();
    before_rebuild(); 
    movingMesh=0;
    conveyor=0;
    skinSafetyFactor=0.;

    //bitmask settings
    EDGE_INACTIVE = new int[3];
    CORNER_INACTIVE = new int[3];
      EDGE_INACTIVE[0] = 1;    EDGE_INACTIVE[1] = 2;    EDGE_INACTIVE[2] = 4;
    CORNER_INACTIVE[0] = 8;  CORNER_INACTIVE[1] = 16; CORNER_INACTIVE[2] = 32;
}

/* ---------------------------------------------------------------------- */

STLtri::~STLtri()
{
   lmp->memory->destroy_3d_double_array(node);
   lmp->memory->destroy_2d_double_array(facenormal);
   lmp->memory->destroy_2d_double_array(f_tri);
   lmp->memory->destroy_2d_double_array(fn_fshear);
   lmp->memory->destroy_2d_double_array(cK);
   lmp->memory->destroy_3d_double_array(ogK);
   lmp->memory->destroy_2d_double_array(ogKlen);
   lmp->memory->destroy_3d_double_array(oKO);
   delete[] rK;
   lmp->memory->destroy_2d_double_array(v);
   lmp->memory->destroy_2d_double_array(xoriginal);
   lmp->memory->destroy_2d_double_array(f);
   delete[] rmass;
   delete[] Area;
   lmp->memory->destroy_2d_int_array(neighfaces);
   delete[] contactInactive;
   
   if(conveyor) lmp->memory->destroy_3d_double_array(v_node);

   delete []EDGE_INACTIVE;delete []CORNER_INACTIVE;
}

/* ---------------------------------------------------------------------- */

void STLtri::grow_arrays()
{
    nTriMax+=DELTATRI;
    xvf_lenMax+=DELTATRI*VECPERTRI;

    node=(double***)(lmp->memory->grow_3d_double_array(node,nTriMax, 3 , 3, "stl_tri_node"));
    node_lastRe=(double***)(lmp->memory->grow_3d_double_array(node_lastRe,nTriMax, 3 , 3, "stl_tri_node_lastRe"));
    facenormal=(double**)(lmp->memory->grow_2d_double_array(facenormal,nTriMax, 3, "stl_tri_facenormal"));
    
    cK=(double**)(lmp->memory->grow_2d_double_array(cK, nTriMax, 3, "stl_tri_cK"));
    ogK=(double***)(lmp->memory->grow_3d_double_array(ogK, nTriMax, 3, 3, "stl_tri_ogK"));
    ogKlen=(double**)(lmp->memory->grow_2d_double_array(ogKlen, nTriMax, 3, "stl_tri_ogKlen"));
    oKO=(double***)(lmp->memory->grow_3d_double_array(oKO, nTriMax, 3, 3, "stl_tri_oKO"));
    rK=(double*)(lmp->memory->srealloc(rK, nTriMax*sizeof(double), "stl_tri_rK"));
    Area=(double*)(lmp->memory->srealloc(Area, nTriMax*sizeof(double), "stl_tri_Area"));
    f_tri=(double**)(lmp->memory->grow_2d_double_array(f_tri,nTriMax,3, "stl_tri_f_tri"));
    fn_fshear=(double**)(lmp->memory->grow_2d_double_array(fn_fshear,nTriMax,3,"stl_tri_fn_fshear"));
    neighfaces=(int**)(lmp->memory->grow_2d_int_array(neighfaces,nTriMax, 3, "stl_tri_neighfaces"));
    contactInactive=(int*)(lmp->memory->srealloc(contactInactive, nTriMax*sizeof(int), "stl_tri_contactActive"));

    for (int i=0;i<nTriMax;i++) {
        f_tri[i][0]=0.;f_tri[i][2]=0.;f_tri[i][2]=0.;
        fn_fshear[i][0]=0.;fn_fshear[i][1]=0.;fn_fshear[i][2]=0.;

        neighfaces[i][0]=-1;neighfaces[i][1]=-1;neighfaces[i][2]=-1;
        contactInactive[i]=0; //init with all contacts active
    }
}

/* ---------------------------------------------------------------------- */
//initialize conveyor model
void STLtri::initConveyor()
{
    conveyor=1;
    if(v_node==NULL) v_node=(double***)(lmp->memory->grow_3d_double_array(v_node,nTriMax, 3 , 3, "stl_tri_v_node"));
    double vtri[3];
    double tmp[3];
    double scp;
    double conv_vel_mag=vectorMag3D(conv_vel);

    for (int i=0;i<nTri;i++)
    {

        scp=vectorDot3D(conv_vel,facenormal[i]);
        vectorScalarMult3D(facenormal[i],scp,tmp);
        for(int j=0;j<3;j++)
        {
            
            vectorSubtract3D(conv_vel,tmp,v_node[i][j]);
            if(vectorMag3D(v_node[i][j])>0.)
            {
                vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
                vectorScalarMult3D(v_node[i][j],conv_vel_mag);
            }
           
        }
    }
}

/* ---------------------------------------------------------------------- */

//copy values of x to xoriginal
void STLtri::initMove(double skinSafety)
{
    if(conveyor) error->all("Can not use moving mesh feature together with conveyor feature");
    movingMesh=1;
    skinSafetyFactor=skinSafety;

    xvf_len=nTri*VECPERTRI;

    if(v_node==NULL) v_node=(double***)(lmp->memory->grow_3d_double_array(v_node,nTriMax, 3 , 3, "stl_tri_v_node"));
    if(x==NULL) x = (double**)(lmp->memory->srealloc(x,xvf_lenMax*sizeof(double *),"stl_tri:x"));
    if(v==NULL) v = (double**)(lmp->memory->grow_2d_double_array(v, xvf_lenMax, 3, "stl_tri:v"));

    //zero node velocity
    for (int i=0;i<xvf_lenMax;i++){
        for (int j=0;j<3;j++) v[i][j]=0.;
    }

    int m=0;
    for (int i=0;i<nTriMax;i++)
    {
        
        for (int j=0;j<3;j++) {
            v_node[i][j]=v[m];
            x[m++]=node[i][j];
        }

        x[m++]=facenormal[i];
        x[m++]=cK[i];
        for (int j=0;j<3;j++) x[m++]=ogK[i][j];
        for (int j=0;j<3;j++) x[m++]=oKO[i][j];
    }
    if(xoriginal==NULL) xoriginal = (double**)(lmp->memory->grow_2d_double_array(xoriginal, xvf_lenMax, 3, "stl_tri:xoriginal"));
    if(f==NULL) f = (double**)(lmp->memory->grow_2d_double_array(f, xvf_lenMax, 3, "stl_tri:f"));
    if(rmass==NULL) rmass=(double*)(lmp->memory->srealloc(rmass, xvf_lenMax*sizeof(double), "stl_import_rmass"));

    vecToPoint();
    for (int i=0;i<xvf_len;i++)
    {
       for(int j=0;j<3;j++) xoriginal[i][j]=x[i][j];
       rmass[i]=1.;
    }
    pointToVec();

}

/* ---------------------------------------------------------------------- */

//save node's positions
void STLtri::before_rebuild()
{
    for (int i=0;i<nTri;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                node_lastRe[i][j][k]=node[i][j][k];
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void STLtri::vecToPoint()
{
    for (int i=0;i<nTri;i++)
    {
        
        for(int j=0;j<3;j++) {
            facenormal[i][j]+=cK[i][j];
        }

        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                oKO[i][j][k]+=node[i][j][k];
                ogK[i][j][k]+=cK[i][k];
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void STLtri::pointToVec()
{
    for (int i=0;i<nTri;i++)
    {
        
        for(int j=0;j<3;j++) facenormal[i][j]-=cK[i][j];

        vectorScalarDiv3D(facenormal[i],vectorMag3D(facenormal[i]));

        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                oKO[i][j][k]-=node[i][j][k];
                ogK[i][j][k]-=cK[i][k];
            }
            
            vectorScalarDiv3D(oKO[i][j],vectorMag3D(oKO[i][j]));
        }
    }

}

/* ---------------------------------------------------------------------- */

void STLtri::pack_restart()
{
   lmp->error->all("Restart not yet enabled for STL_tri class");
}

/* ---------------------------------------------------------------------- */

void STLtri::unpack_restart()
{
   lmp->error->all("Restart not yet enabled for STL_tri class");
}

/* ---------------------------------------------------------------------- */

void STLtri::getTriBoundBox(int iTri, double *xmin, double *xmax,double delta)
{
    
    for(int j=0;j<3;j++)
    {
      xmin[j]=fmin(node[iTri][0][j],fmin(node[iTri][1][j],node[iTri][2][j]))-delta;
      xmax[j]=fmax(node[iTri][0][j],fmax(node[iTri][1][j],node[iTri][2][j]))+delta;
    }
}

/* ---------------------------------------------------------------------- */

bool STLtri::are_coplanar_neighs(int i1,int i2)
{
    
    if     (neighfaces[i1][0] == i2)
        return contactInactive[i1] & EDGE_INACTIVE[0];
    else if(neighfaces[i1][1] == i2)
        return contactInactive[i1] & EDGE_INACTIVE[1];
    else if(neighfaces[i1][2] == i2)
        return contactInactive[i1] & EDGE_INACTIVE[2];
    else return false;
}
