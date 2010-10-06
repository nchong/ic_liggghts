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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_meshGran.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "triSpherePenetration.h"
#include "stl_tri.h"
#include "mpi.h"

using namespace LAMMPS_NS;

#define EPSILON 0.0001
#define SIDE_RATIO_TOLERANCE 10
#define NEIGH_TOLERANCE 1e-5

/* ---------------------------------------------------------------------- */

FixMeshGran::FixMeshGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
   if ((!atom->radius_flag)||(!atom->rmass_flag)) lmp->error->all("Particles need to store radius and mass for granular tri-particle interaction, use atom style granular");
  if (domain->box_exist==0) error->all("Must define simulation box before using 'fix mesh/gran'");

  if (narg < 12) error->all("Illegal fix mesh/gran command, not enough arguments");

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all("Fix mesh/gran requires atom attributes radius, rmass (atom style granular)");

  stl_filename = arg[3];
  atom_type_wall = force->inumeric(arg[4]);
  scale_fact=force->numeric(arg[5]);

  off_fact=new double[3];
  off_fact[0]=force->numeric(arg[6]);
  off_fact[1]=force->numeric(arg[7]);
  off_fact[2]=force->numeric(arg[8]);

  rot_angle=new double[3];
  rot_angle[0]=force->numeric(arg[9]);
  rot_angle[1]=force->numeric(arg[10]);
  rot_angle[2]=force->numeric(arg[11]);

  iarg=12;

  nevery=1; 

  if (comm->me == 0) {
    fprintf(screen,"\nImporting STL file '%s' \n",stl_filename);
  }

  STLdata= new STLtri(lmp);

  //allocate input class
  mystl_input = new Input(lmp,0,NULL);

  mystl_input->stlfile(stl_filename,this);

  nTri=STLdata->nTri;
  node=STLdata->node;
  EDGE_INACTIVE=STLdata->EDGE_INACTIVE;
  CORNER_INACTIVE=STLdata->CORNER_INACTIVE;

  calcTriCharacteristics(STLdata->nTri,STLdata->node,STLdata->cK,STLdata->ogK,STLdata->ogKlen,STLdata->oKO,STLdata->rK,STLdata->Area,STLdata->facenormal,STLdata->neighfaces,STLdata->contactInactive);

  if  (comm->me==0) fprintf(screen,"\nImport of %d triangles completed successfully!\n\n",STLdata->nTri);

  force_total=new double[3];
  torque_total=new double[3];
  p_ref=new double[3];

  analyseStress=false;

  force_total[0]=0.;force_total[1]=0.;force_total[2]=0.;
  torque_total[0]=0.;torque_total[2]=0.;torque_total[2]=0.;
  p_ref[0]=0.;p_ref[2]=0.;p_ref[2]=0.;

  restart_global=1;

  STLdata->conv_vel[0]=0.;
  STLdata->conv_vel[1]=0.;
  STLdata->conv_vel[2]=0.;

  if(narg>=16 && strcmp(arg[iarg],"conveyor")==0)
  {
      iarg++;
      STLdata->conv_vel[0]=atof(arg[iarg++]);
      STLdata->conv_vel[1]=atof(arg[iarg++]);
      STLdata->conv_vel[2]=atof(arg[iarg++]);
      double conv_vel_mag=vectorMag3D(STLdata->conv_vel);
      if(conv_vel_mag<EPSILON) error->all("Conveyor velocity too low");
      STLdata->initConveyor();
  }
}

/* ---------------------------------------------------------------------- */

void FixMeshGran::calcTriCharacteristics(int nTri,double ***node,double **cK,double ***ogK,double **ogKlen,double ***oKO,double *rK,double *Area,double**facenormal,int **neighfaces,int *contactInactive)
{
      double temp[3];
      bool flag_skewed = false;

      p_ref = new double[3];
      p_ref[0]=0.;p_ref[1]=0.;p_ref[2]=0.;
      double A_tri,A_ges=0.;
      double oKO0Neg[3];
      int count = 0;

      for(int i=0;i<nTri;i++)
      {
          for (int j=0;j<3;j++){  //loop components
            oKO[i][0][j]=node[i][1][j]-node[i][0][j];
            oKO[i][1][j]=node[i][2][j]-node[i][1][j];
            oKO[i][2][j]=node[i][0][j]-node[i][2][j];
          }
          vectorScalarMult3D(oKO[i][0],-1.,oKO0Neg);
          vectorCross3D(oKO0Neg,oKO[i][1],temp);
          A_tri=0.5*vectorMag3D(temp);
          Area[i]=A_tri;

          vectorCopy3D(node[i][0],temp);
          vectorAdd3D(temp,node[i][1],temp);
          vectorAdd3D(temp,node[i][2],temp);
          vectorScalarDiv3D(temp,3.); 

          vectorScalarMult3D(p_ref,A_ges);
          vectorScalarMult3D(temp,A_tri);
          vectorAdd3D(p_ref,temp,p_ref);
          A_ges+=A_tri;
          vectorScalarDiv3D(p_ref,A_ges);
      }
      
      for(int i=0;i<nTri;i++) //loop tris
      {
          for (int j=0;j<3;j++){  //loop components
            oKO[i][0][j]=node[i][1][j]-node[i][0][j];
            oKO[i][1][j]=node[i][2][j]-node[i][1][j];
            oKO[i][2][j]=node[i][0][j]-node[i][2][j];
          }

          double mag0=vectorMag3D(oKO[i][0]),mag1=vectorMag3D(oKO[i][1]),mag2=vectorMag3D(oKO[i][2]);
          double sr0 = fmax(mag0/mag1,mag1/mag0),sr1 = fmax(mag1/mag2,mag2/mag1),sr2 = fmax(mag0/mag2,mag2/mag0);
          double side_ratio_max = fmax(fmax(sr0,sr1),sr2);
          if(side_ratio_max > SIDE_RATIO_TOLERANCE) flag_skewed=true;

          //normalize
          for (int j=0;j<3;j++) vectorScalarDiv3D(oKO[i][j],vectorMag3D(oKO[i][j]));

          double *ans=new double[3], *v=new double[3], *add_v=new double[3];
          v[0]=0.;v[1]=0.;v[2]=0.;ans[0]=0.;ans[1]=0.;ans[2]=0.;
          add_v[0]=0.;add_v[1]=0.;add_v[2]=0.;

          double M[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
          double Mt[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

          vectorSubtract3D(node[i][0],node[i][2],M[0]);
          vectorSubtract3D(oKO[i][0],oKO[i][2],M[1]);
          vectorSubtract3D(oKO[i][1],oKO[i][2],M[2]);

          vectorCross3D(M[1],M[2],add_v);
          vectorAdd3D(M[0],add_v,M[0]);
          vectorAdd3D(v,add_v,v);

          //transpose matrix
          MathExtra::transpose3(M,Mt);

          MathExtra::mldivide3(Mt, v, ans, lmp->error);

          //ans[0] should be equal to 1, otherwise the operation has failed
          if((fabs(ans[0]-1)>1e-3)||std::isnan(ans[0]))
          {

              char *buff=new char[100];
              sprintf(buff,"STL import failed at triangle #%d (line %d), degenerated triangle?",i,i*7+1);
              lmp->error->all(buff);
          }

          double lambda = ans[1]/ans[0];
          double mu = ans[2]/ans[0];

          vectorSubtract3D(oKO[i][0],oKO[i][2],cK[i]);
          vectorScalarMult3D(cK[i],lambda);
          vectorAdd3D(cK[i],node[i][0],cK[i]);

          for(int j=0;j<3;j++)
          {
              vectorSubtract3D(cK[i],node[i][j],temp);
              double dot=vectorDot3D(temp,oKO[i][j]);
              vectorScalarMult3D(oKO[i][j],dot,temp);
              vectorCopy3D(temp,ogK[i][j]);

              vectorAdd3D(ogK[i][j],node[i][j],ogK[i][j]);
              vectorSubtract3D(ogK[i][j],cK[i],ogK[i][j]);

              ogKlen[i][j]=vectorMag3D(ogK[i][j]);
          }

          rK[i]=vectorMag3D(ogK[i][0]);

      }
      if(flag_skewed && comm->me ==0) error->warning("Imported mesh contains skewed triangles with a size ratio larger 10");

      if(comm->me==0)fprintf(screen,"Mesh calculations running. This may take a while...");
      int common_c, common1[3], common2[3];
      for(int i=0;i<nTri;i++) 
      {
          for(int j=i+1;j<nTri;j++) 
          {
              
              common_c=0;
              for(int k=0;k<3;k++)
              {
                  common1[k]=0;common2[k]=0;
              }
              
              for(int iVertex=0;iVertex<3;iVertex++)
              {
                 for(int jVertex=0;jVertex<3;jVertex++)
                 {
                     vectorSubtract3D(node[i][iVertex],node[j][jVertex],temp);
                     
                     double dist=vectorMag3D(temp);
                     if(dist < NEIGH_TOLERANCE*0.5*(rK[i]+rK[j]))
                     {
                         
                         common1[iVertex]=1;
                         common2[jVertex]=1;
                         common_c++;
                     }
                 }
              }
              
              if(common_c==3)
              {
                  char *buff=new char[100];
                  sprintf(buff,"STL import failed because triangles #%d (line %d) and #%d (line %d) are duplicate",i,i*7+1,j,j*7+1);
                  lmp->error->all(buff);
              }
              
              if(common_c==2)
              {
                  
                  if(common1[0]&&common1[1]) neighfaces[i][0] = j;
                  if(common1[1]&&common1[2]) neighfaces[i][1] = j;
                  if(common1[2]&&common1[0]) neighfaces[i][2] = j;

                  if(common2[0]&&common2[1]) neighfaces[j][0] = i;
                  if(common2[1]&&common2[2]) neighfaces[j][1] = i;
                  if(common2[2]&&common2[0]) neighfaces[j][2] = i;
              }
          }
      }

      int *prevFaces = new int[nTri];

      for(int i = 0; i < nTri; i++) //loop tris
      {
          
          for(int j=0;j<3;j++)
          {
              int iNeigh=neighfaces[i][j];
              if(iNeigh >= 0)
              {
                  
                  double dot = vectorDot3D(facenormal[i],facenormal[iNeigh]);
                  if(fabs(dot) > (1.-EPSILON) )
                  {
                      contactInactive[i] |= EDGE_INACTIVE[j];
                      
                  }

                  if(iNeigh > i)
                  {
                      contactInactive[i] |= EDGE_INACTIVE[j];
                      
                  }
              }
          }
          
          for(int j=0;j<3;j++)
          {
              int edge1 = j;
              int edge2 = j-1;
              if(edge2 < 0) edge2=2;

              int neigh1 = neighfaces[i][edge1];
              int neigh2 = neighfaces[i][edge2];

              if(neigh1 >= 0 && neigh2 >= 0)
              {
                  double dot1=vectorDot3D(facenormal[i],facenormal[neigh1]);
                  double dot2=vectorDot3D(facenormal[i],facenormal[neigh2]);
                  if(fabs(dot1)>(1-EPSILON) && fabs(dot2)>(1-EPSILON))
                  {
                      contactInactive[i] |= CORNER_INACTIVE[j];
                      
                  }
              }

              int maxid = -1,nPrev = 0;
              
              if(neigh1 >= 0) maxid = get_max_index_sharedcorner(i,nPrev,prevFaces,node[i][j],node,rK,neighfaces);

              if(i < maxid)
              {
                  contactInactive[i] |= CORNER_INACTIVE[j];
                  
              }
          }
      }
      delete []prevFaces;
      
      if(comm->me==0)fprintf(screen,"finished!\n");
}

/* ---------------------------------------------------------------------- */

int FixMeshGran::get_max_index_sharedcorner(int iTri,int &nPrev,int *prevFaces,double *node2check,double ***node,double *rK,int **neighfaces)
{
    
    if(iTri < 0)
    {
        
        return -1;
    }

    double temp[3];

    for(int ic = 0; ic < nPrev-1; ic++)
       if(prevFaces[ic] == iTri) return -1;

    prevFaces[nPrev++] = iTri;

    bool contains = false;
    for(int j = 0; j < 3; j++)
    {
        vectorSubtract3D(node2check,node[iTri][j],temp);
        double dist=vectorMag3D(temp);
        if(dist < NEIGH_TOLERANCE*rK[iTri]) contains = true;
    }

    if(!contains)
    {
        
        return -1;
    }

    int n0 = neighfaces[iTri][0], n1 = neighfaces[iTri][1], n2 = neighfaces[iTri][2];

    int n0max = -1, n1max = -1, n2max = -1;

    n0max = get_max_index_sharedcorner(n0,nPrev,prevFaces,node2check,node,rK,neighfaces);
    n1max = get_max_index_sharedcorner(n1,nPrev,prevFaces,node2check,node,rK,neighfaces);
    n2max = get_max_index_sharedcorner(n2,nPrev,prevFaces,node2check,node,rK,neighfaces);

    //return the maximum index including this tri and its 3 neighs
    if(iTri > n0max && iTri > n1max && iTri > n2max) return iTri;
    else if(n0max > n1max && n0max > n2max) return n0max;
    else if(n1max > n2max ) return n1max;
    else return n2max;
}

/* ---------------------------------------------------------------------- */

FixMeshGran::~FixMeshGran()
{
   delete mystl_input;
   delete STLdata;
   delete []force_total;
   delete []torque_total;
   delete []p_ref;
}

/* ---------------------------------------------------------------------- */

int FixMeshGran::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMeshGran::write_restart(FILE *fp)
{
  int n = 0;

  int listlen=12*nTri+2;
  double *list=new double[listlen];

  list[n++]=static_cast<double>(STLdata->nTri);

  for(int i=0;i<nTri;i++)
  {
    for(int j=0;j<3;j++) list[n++] = STLdata->facenormal[i][j];

    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        list[n++] = STLdata->node[i][j][k];
  }

  if (comm->me == 0) {
    int size = (n+n_children()) * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
    children_write(fp);
  }
  delete []list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMeshGran::restart(char *buf)
{

  int n = 0;
  double *list = (double *) buf;

  int nTri_restart = static_cast<int> (list[n++]);

  while(STLdata->nTriMax<nTri_restart) STLdata->grow_arrays();
  STLdata->nTri=nTri_restart;
  STLdata->xvf_len=STLdata->nTri*VECPERTRI;

  for(int i=0;i<nTri;i++)
  {
    for(int j=0;j<3;j++) STLdata->facenormal[i][j] = list[n++];

    for(int j=0;j<3;j++)
      for(int k=0;k<3;k++)
        STLdata->node[i][j][k] = list[n++];
  }

  nTri=STLdata->nTri;
  node=STLdata->node;
  calcTriCharacteristics(STLdata->nTri,STLdata->node,STLdata->cK,STLdata->ogK,STLdata->ogKlen,STLdata->oKO,STLdata->rK,STLdata->Area,STLdata->facenormal,STLdata->neighfaces,STLdata->contactInactive);

  children_restart(&(list[n]));
}
