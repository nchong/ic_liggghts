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
/* NP ----------------------------------------------------------------------
   The algorithms are partly based on the follwing papers:
   1: Int. J. Numer. Meth. Engng 2001, 51, 1407-1421
   2: Int. J. Numer. Meth. Engng 2001, 51, 1423-1436
   equation numbers refer to paper 1
------------------------------------------------------------------------- */
#ifndef LMP_TRISPHEREPENETRATION_H
#define LMP_TRISPHEREPENETRATION_H

#include "lammps.h"
#include "fix.h"
#include "input.h"
#include "math.h"
#include "math_extra.h"
#include "myvector.h"
#include "stl_tri.h"
#include "atom.h"
#include "lammps.h"
#include "error.h"
#include "update.h"

#define LARGE 10.

using namespace LAMMPS_NS;

namespace TRISPHERE {
    double resolveTriSphereContact(class LAMMPS *lmp,int iPart,int iTri,STLtri* STL_tri,double fSh,double *enO,bool neighbuild,double treshold,bool &excl);
    double returnOverlap(class LAMMPS *lmp,double *cs,double *xc,double *t,double rs,double *enO,bool neighbuild);
};

  inline double TRISPHERE::resolveTriSphereContact(class LAMMPS *lmp,int iPart,int iTri,STLtri* STL_tri,double fSh,double *enO,bool neighbuild,double treshold,bool &excl)
  {
      //returns value <0 for contact, and value >0 for no contact
      //also returns value for enO (normalized vector from the spheres center to the contact point with the tri)
      //if the "neighbuild flag" is set [for neigh list build], the distance is returned - >0 for no overlap and in case of overlap, the overklap <0
      //if the "neighbuild" flag is not set [for force calculation], a value >0 (but not equal to the distance) is returned for no contact, and the overlap <0 for the overlap case
      
      //allocate memory
      double xc[3];
      double vxc[3];
      double t[3];
      double t2[3];
      double bi[3];
      double pxcPrev[3];
      double pxcNext[3];

      double *ck,*noKO;  //center of tri, face normal
      double **ogK,*ogKlen,**nodes;
      double rs,*cs;   //sphere (particle) radius, sphere center
      double l; //penetration

      ck=STL_tri->cK[iTri];
      noKO=STL_tri->facenormal[iTri];

      rs=lmp->atom->radius[iPart];
      cs=lmp->atom->x[iPart];

      vectorSubtract3D(cs,ck,t);
      double a = vectorDot3D(t,noKO);
      l=fabs(a)-rs;

      bool surfacePen=true;
      excl=false;

      if (l>0)
      {
          
          surfacePen=false;
          if (!neighbuild)  return l;
          else if(l>treshold) return l;
      }
      
      nodes=STL_tri->node[iTri];
      ogK=STL_tri->ogK[iTri];
      ogKlen=STL_tri->ogKlen[iTri];
      int contactInactive=STL_tri->contactInactive[iTri];
      int *EDGE_INACTIVE=STL_tri->EDGE_INACTIVE;
      int *CORNER_INACTIVE=STL_tri->CORNER_INACTIVE;

      //start checking
      vectorCopy3D(noKO,xc);
      vectorScalarMult3D(xc,a);
      vectorSubtract3D(cs,xc,xc);
      vectorSubtract3D(xc,ck,vxc);

      double vxcmag=vectorMag3D(vxc);
      for (int i=0;i<3;i++) bi[i]=vectorDot3D(vxc,ogK[i])/(vxcmag*ogKlen[i]); //equ 18

      int iSub;
      
      if (bi[0]>bi[1] && bi[0]>bi[2]) iSub=0;
      else if (bi[1]>bi[2]) iSub=1;
      else iSub=2;

      if (bi[iSub]*bi[iSub]*vxcmag*vxcmag+2*fSh*rs*bi[iSub]*vxcmag+fSh*rs*fSh*rs<(vectorMag3DSquared(ogK[iSub])))
      {
          
          return returnOverlap(lmp,cs,xc,t,rs,enO,neighbuild);
      }

      if (bi[iSub]*vxcmag>(ogKlen[iSub]+rs))
      {
          
          if(neighbuild) //calc new xc, vxc to be able to return exact distance to tri
          {
             vectorSubtract3D(vxc,ogK[iSub],t);
             double scalProd=vectorDot3D(t,ogK[iSub])/ogKlen[iSub];
             vectorScalarMult3D(ogK[iSub],scalProd/ogKlen[iSub],t);
             vectorSubtract3D(xc,t,xc);
             vectorSubtract3D(xc,ck,vxc);  
             return returnOverlap(lmp,cs,xc,t,rs,enO,neighbuild);
             
          }
          else return LARGE;

      }

      if(!neighbuild && contactInactive==63)
      {
          
          return LARGE;
      }
      else if (contactInactive==63) excl=true;

      vectorSubtract3D(vxc,ogK[iSub],t);
      double scalProd=vectorDot3D(t,ogK[iSub])/ogKlen[iSub];
      vectorScalarMult3D(ogK[iSub],scalProd/ogKlen[iSub],t);
      vectorSubtract3D(xc,t,xc);
      vectorSubtract3D(xc,ck,vxc);  
      
      //get values for i-1 and i+1

      int iSubPrev=iSub-1;
      if (iSubPrev<0) iSubPrev=2;
      int iSubNext=iSub+1;
      if (iSubNext>2) iSubNext=0;

      vectorScalarDiv3D(ogK[iSubPrev],ogKlen[iSubPrev],pxcPrev);
      vectorScalarMult3D(pxcPrev,vectorDot3D(vxc,pxcPrev),pxcPrev);         

      vectorScalarDiv3D(ogK[iSubNext],ogKlen[iSubNext],pxcNext);
      vectorScalarMult3D(pxcNext,vectorDot3D(vxc,pxcNext),pxcNext);         

      vectorSubtract3D(cs,nodes[iSub],t);
      vectorSubtract3D(cs,nodes[iSubNext],t2);

      bool cornerI     = ((vectorMag3DSquared(pxcPrev)>vectorMag3DSquared(ogK[iSubPrev]))&&(vectorMag3DSquared(t )<rs*rs));
      bool cornerNextI = ((vectorMag3DSquared(pxcNext)>vectorMag3DSquared(ogK[iSubNext]))&&(vectorMag3DSquared(t2)<rs*rs));

      if (cornerI || cornerNextI )  
      {
          
          if (cornerI) vectorCopy3D(nodes[iSub],xc);
          else if (cornerNextI) vectorCopy3D(nodes[iSubNext],xc);

          if(( cornerI && contactInactive&CORNER_INACTIVE[iSub]) || (cornerNextI && contactInactive&CORNER_INACTIVE[iSubNext]))
          {
              
              if (!neighbuild )return LARGE;
              excl = true;
          }
          return returnOverlap(lmp,cs,xc,t,rs,enO,neighbuild);
      }

      if ((vectorMag3D(t)-rs>=0.) || (vectorMag3D(t2)-rs>=0.))
      {
          
      }

      if ((vectorMag3DSquared(pxcPrev)<=vectorMag3DSquared(ogK[iSubPrev])||(vectorDot3D(pxcPrev,ogK[iSubPrev])<0)) && (vectorMag3DSquared(pxcNext)<=vectorMag3DSquared(ogK[iSubNext])||(vectorDot3D(pxcNext,ogK[iSubNext])<0)))
      {
          
          if(contactInactive&EDGE_INACTIVE[iSub])
          {
              
              if(!neighbuild) return LARGE;
              excl = true;
          }
          return returnOverlap(lmp,cs,xc,t,rs,enO,neighbuild);
      }
      else 
      {
          
          if (neighbuild)
          {
              //set xc as corner with the minimum distance
              vectorSubtract3D(nodes[iSub],cs,t);
              vectorSubtract3D(nodes[iSubNext],cs,t2);
              if (vectorMag3D(t)<vectorMag3D(t2)) vectorCopy3D(nodes[iSub],xc);
              else vectorCopy3D(nodes[iSubNext],xc);
              return returnOverlap(lmp,cs,xc,t,rs,enO,neighbuild);
          }
          else return LARGE; //just return 1 to signal "no contact"
      }
  }

  inline double TRISPHERE::returnOverlap(class LAMMPS *lmp,double *cs,double *xc,double *t,double rs,double *enO,bool neighbuild)
  {
      vectorSubtract3D(cs,xc,t);
      if(!neighbuild)
      {
          vectorSubtract3D(xc,cs,enO); //equ. 32
          vectorScalarDiv3D(enO,vectorMag3D(enO),enO);
      }
      double overl=(vectorMag3D(t)-rs);//equ. 31

      return overl;
  }
#endif
