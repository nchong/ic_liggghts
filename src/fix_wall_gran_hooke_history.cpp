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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran_hooke_history.h"
#include "pair_gran_hooke_history.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "triSpherePenetration.h"
#include "fix_propertyGlobal.h"
#include "fix_meshGran.h"
#include "fix_tri_neighlist.h"
#include "comm.h"
#include "mech_param_gran.h"
#include "fix_rigid.h"
#include "debug_single_CAD.h"

using namespace LAMMPS_NS;

enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

#define DELTA_TRI_CONTACTS 6 

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define EPSILON_MOVINGMESH 1e-12

/* ---------------------------------------------------------------------- */

FixWallGranHookeHistory::FixWallGranHookeHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if(strncmp(style,"wall/gran/hooke",15)!=0  && strncmp(style,"wall/gran/hertz",15)!=0) return;

  if (narg < 7) error->all("Illegal fix wall/gran command");

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Fix wall/gran requires atom attributes radius, omega, torque");

  FixMeshGranList=NULL;
  fix_tri_neighlist=NULL;

  // wall/particle coefficients
  dampflag = atoi(arg[3]);
  cohesionflag = atoi(arg[4]);

  lo = hi = 0.0;

  if (cohesionflag < 0 || cohesionflag > 1 || dampflag < 0 || dampflag > 1)
    error->all("Illegal fix wall/gran command");

  int iarg = 5;
  if (strncmp(arg[iarg],"mesh/gran",8) == 0) {
      wallstyle=MESHGRAN_FWGHH;

      for (int ifix=0;ifix<modify->nfix;ifix++)
      {
          if(strncmp(modify->fix[ifix]->style,"wall/gran/",9)==0)
          {
              FixWallGranHookeHistory* othF=static_cast<FixWallGranHookeHistory*>(modify->fix[ifix]);
              if (othF->wallstyle==MESHGRAN_FWGHH) error->all("There must be only one fix wall/gran with style 'mesh/gran'");
          }
      }

      if (narg < iarg+2) error->all("Illegal fix wall/gran command, not enough arguments");
      nFixMeshGran = atoi(arg[iarg+1]);
      if (narg < iarg+2+nFixMeshGran) error->all("Illegal fix wall/gran command, not enough arguments");
      FixMeshGranList=new FixMeshGran*[nFixMeshGran];
      for (int i=0;i<nFixMeshGran;i++) FixMeshGranList[i]=static_cast<FixMeshGran*>(NULL);
      for(int i=0;i<nFixMeshGran;i++)
      {
          int f_i=modify->find_fix(arg[iarg+2+i]);
          if (f_i==-1) error->all("Could not find fix mesh/gran id you provided for the fix wall/gran command");
          if (strncmp(modify->fix[f_i]->style,"mesh/gran",8) != 0) error->all("fix wall/gran: A fix belonging to the id you provided is not of type mesh/gran");
          FixMeshGranList[i]=static_cast<FixMeshGran*>(modify->fix[f_i]);
      }

      iarg += (2+nFixMeshGran);
  }else if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+4) error->all("Illegal fix wall/gran command");
    wallstyle = XPLANE_FWGHH;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    atom_type_wall=atoi(arg[iarg+3]);
    iarg += 4;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/gran command");
    wallstyle = YPLANE_FWGHH;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    atom_type_wall=atoi(arg[iarg+3]);
    iarg += 4;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/gran command");
    wallstyle = ZPLANE_FWGHH;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = atof(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = atof(arg[iarg+2]);
    atom_type_wall=atoi(arg[iarg+3]);
    iarg += 4;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/gran command");
    wallstyle = ZCYLINDER_FWGHH;
    lo = hi = 0.0;
    cylradius = atof(arg[iarg+1]);
    atom_type_wall=atoi(arg[iarg+2]);
    iarg += 3;
  }

  if(lo > hi) error->all("Fix wall/gran: lo value of wall position must be < hi value");

  // check for trailing keyword/values

  wiggle = 0;
  wshear = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"wiggle") == 0) {
      if (wallstyle==MESHGRAN_FWGHH) error->all("Fix wall/gran: Can not use wiggle together with style mesh/gran. Please use fix move/mesh/gran for moving mesh");
      if (iarg+4 > narg) error->all("Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all("Illegal fix wall/gran command");
      amplitude = atof(arg[iarg+2]);
      period = atof(arg[iarg+3]);
      wiggle = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"shear") == 0) {
      if (wallstyle==MESHGRAN_FWGHH) error->all("Fix wall/gran: Can not use shear together with style mesh/gran. Please use fix move/mesh/gran for moving mesh");
      if (iarg+3 > narg) error->all("Illegal fix wall/gran command");
      if (strcmp(arg[iarg+1],"x") == 0) axis = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) axis = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) axis = 2;
      else error->all("Illegal fix wall/gran command");
      vshear = atof(arg[iarg+2]);
      wshear = 1;
      iarg += 3;
    } else error->all("Illegal fix wall/gran command");
  }

  if (wallstyle == XPLANE_FWGHH && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE_FWGHH && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE_FWGHH && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER_FWGHH && (domain->xperiodic || domain->yperiodic))
    error->all("Cannot use wall in periodic dimension");

  if (wiggle && wshear) error->all("Cannot wiggle and shear fix wall/gran");
  if (wiggle && wallstyle == ZCYLINDER_FWGHH && axis != 2)
    error->all("Invalid wiggle direction for fix wall/gran");
  if (wshear && wallstyle == XPLANE_FWGHH && axis == 0)
    error->all("Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == YPLANE_FWGHH && axis == 1)
    error->all("Invalid shear direction for fix wall/gran");
  if (wshear && wallstyle == ZPLANE_FWGHH && axis == 2)
    error->all("Invalid shear direction for fix wall/gran");

  // setup oscillations

  if (wiggle) {
    double PI = 4.0 * atan(1.0);
    omega = 2.0*PI / period;
  }

  restart_global = 1;
  restart_peratom = 1;
  create_attribute = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  shearhistory=1;
  shear = NULL;
  npartners = NULL;
  partner = NULL;

  if (wallstyle==MESHGRAN_FWGHH) maxpartners=DELTA_TRI_CONTACTS;
  else maxpartners=1;

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  initSubstyle();
  time_depend = 1;

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixWallGranHookeHistory::~FixWallGranHookeHistory()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays
  memory->destroy_3d_double_array(shear);
  memory->destroy_3d_int_array(partner);
  delete []npartners;

}
/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::initSubstyle()
{
  
  if (strcmp(this->style,"wall/gran/hooke/history")!=0) return;

  // initialize as if particle is not touching wall
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
  {
      npartners[i]=0;
      for (int k=0;k<maxpartners;k++)
      {
          partner[i][k][0] = partner[i][k][0] = -1;
          shear[i][k][0] = shear[i][k][1] = shear[i][k][2] = 0.0;
      }
  }

  if (!(force->pair_match("gran/hooke/history",0))) error->all("Fix wall/gran/hooke/history can only be used together with pair style gran/hooke/history");

  mpg = static_cast<PairGranHookeHistory*>(force->pair)->getMatProp();
}

/* ---------------------------------------------------------------------- */

int FixWallGranHookeHistory::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::registerTriNeighlist(int miflag)
{
    char **fixarg;
    fixarg=new char*[4];
    for (int kk=0;kk<4;kk++) fixarg[kk]=new char[30];
    //register  neighborlist fix
    fixarg[0]="fix_neigh_tri";
    fixarg[1]="all";
    fixarg[2]="neighlist/tri";
    strcpy(fixarg[3],id);
    modify->add_fix(4,fixarg);
    delete []fixarg;

    int i_fix=modify->find_fix("fix_neigh_tri");
    if (i_fix==-1) error->all("Could not register a fix for the triangle neighbor list");

    fix_tri_neighlist=static_cast<FixTriNeighlist*>(modify->fix[i_fix]);

    if(miflag) modify->init();
    
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::init()
{
  
  dt = update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  //register a fix TriNeighlist
  if ((wallstyle==MESHGRAN_FWGHH) && (fix_tri_neighlist==NULL)) registerTriNeighlist(1);

  if (!(force->pair_match("gran/hooke/history",0))) error->all("Fix wall/gran/hooke/history can only be used together with pair style gran/hooke/history");

  //get material properties
  mpg = static_cast<PairGranHookeHistory*>(force->pair)->getMatProp();

  //check if a fix rigid is registered - important for damp
  fr=NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
  {
      if(strncmp(modify->fix[ifix]->style,"rigid",5)==0) fr=static_cast<FixRigid*>(modify->fix[ifix]);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::setup(int vflag)
{
  
  if ((wallstyle==MESHGRAN_FWGHH) && (fix_tri_neighlist==NULL)) registerTriNeighlist(1);

  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::post_force(int vflag)
{
  double vwall[3],dx,dy,dz,del1,del2,delxy,delr,rsq;
  double kn,kt,gamman,gammat,xmu;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly
  
  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    wlo = lo + amplitude*sin(arg);//amplitude - amplitude*cos(arg);
    whi = hi + amplitude*sin(arg);//amplitude - amplitude*cos(arg);
    vwall[axis] = amplitude*omega*cos(arg);//amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set shear if pair potential stores history

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *shr;

  int iTri,iFMG;
  double deltan;
  double en0[3];
  bool excl;
  int *nTriList;
  int ***tri_neighlist;
  int *delflag;

  int iContactList;

  double contactPoint[3];
  double contactPoint_triCoo[3];
  double **vnode;
  double Mt[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  double wallforce[3];
  double tmp[3];
  double tmp2[3];

  if(wallstyle == MESHGRAN_FWGHH){
      if(fix_tri_neighlist==NULL) error->all("fatal: Could not find triangle neighbor list");
      nTriList=fix_tri_neighlist->nTriList;
      tri_neighlist=fix_tri_neighlist->tri_neighlist;
      delflag=fix_tri_neighlist->delflag;

      //zero wall forces
      reset_wall_forces();
      
  }
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE_FWGHH) {
	    del1 = x[i][0] - wlo;
	    del2 = whi - x[i][0];
	    if (del1 < del2) dx = del1;
	    else dx = -del2;
      } else if (wallstyle == YPLANE_FWGHH) {
	    del1 = x[i][1] - wlo;
	    del2 = whi - x[i][1];
	    if (del1 < del2) dy = del1;
	    else dy = -del2;
      } else if (wallstyle == ZPLANE_FWGHH) {
	    del1 = x[i][2] - wlo;
	    del2 = whi - x[i][2];
	    if (del1 < del2) dz = del1;
	    else dz = -del2;
      } else if (wallstyle == ZCYLINDER_FWGHH) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
	    if (delr > radius[i]) dz = cylradius;
	    else {
	      dx = -delr/delxy * x[i][0];
	      dy = -delr/delxy * x[i][1];
          if (wshear && axis != 2) {
	        vwall[0] = vshear * x[i][1]/delxy;
	        vwall[1] = -vshear * x[i][0]/delxy;
	        vwall[2] = 0.0;
	      }
	     }
      }

      //for style !=MESHGRAN, we have just one wall - handle it as we did in the old LAMMPS wall handling
      if (wallstyle != MESHGRAN_FWGHH){
          rsq = dx*dx + dy*dy + dz*dz;
          if (rsq > radius[i]*radius[i])  resetShearHistory(i,0);
          else
          {
             double deltan=radius[i]-sqrt(rsq);
             deriveContactModelParams(i,deltan,kn,kt,gamman,gammat,xmu);
             
             shr=NULL;
             if (shear!=NULL) shr=shear[i][0];
             
             compute_force(i,rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],radius[i],rmass[i],shr,kn,kt,gamman,gammat,xmu);
            
          }
      }
      else //wallstyle == MESHGRAN
      {
          
          for (int j=0;j<nTriList[i];j++)
          {
              
              delflag[j] = 0;

              iFMG=tri_neighlist[i][j][0];
              iTri=tri_neighlist[i][j][1];

              deltan=TRISPHERE::resolveTriSphereContact(lmp,i,iTri,FixMeshGranList[iFMG]->STLdata,F_SHRINKAGE,en0,false,0.,excl);

              if (deltan<0.) //overlap
              {
                  //add to list of contacts
                  iContactList = add_to_contact_list(i,iFMG,iTri);

                  vectorScalarMult3D(en0,atom->radius[i],contactPoint);
                  vectorAdd3D(atom->x[i],contactPoint,contactPoint);

                  if(FixMeshGranList[iFMG]->STLdata->movingMesh)
                  {
                      
                      for (int m1=0;m1<3;m1++){
                          for (int m2=0;m2<3;m2++){
                              Mt[m1][m2]=FixMeshGranList[iFMG]->STLdata->node[iTri][m2][m1];
                          }
                      }
                      
                      if(fabs(MathExtra::mdet(Mt,lmp->error))>EPSILON_MOVINGMESH) MathExtra::mldivide3(Mt,contactPoint,contactPoint_triCoo,lmp->error);

                      else  
                      {
                          vectorCross3D(FixMeshGranList[iFMG]->STLdata->node[iTri][0],FixMeshGranList[iFMG]->STLdata->node[iTri][1],tmp);
                          if (vectorMag3D(tmp)>EPSILON_MOVINGMESH) for(int m1=0;m1<3;m1++) Mt[m1][2]=tmp[m1];
                          else {
                              vectorCross3D(FixMeshGranList[iFMG]->STLdata->node[iTri][0],FixMeshGranList[iFMG]->STLdata->node[iTri][2],tmp);
                              if (vectorMag3D(tmp)>EPSILON_MOVINGMESH) for(int m1=0;m1<3;m1++) Mt[m1][1]=tmp[m1];
                              else {
                                 vectorCross3D(FixMeshGranList[iFMG]->STLdata->node[iTri][1],FixMeshGranList[iFMG]->STLdata->node[iTri][2],tmp);
                                 if (vectorMag3D(tmp)>EPSILON_MOVINGMESH) for(int m1=0;m1<3;m1++) Mt[m1][0]=tmp[m1];
                              }
                          }
                          MathExtra::mldivide3(Mt,contactPoint,contactPoint_triCoo,lmp->error);
                      }

                      if(std::isnan(contactPoint_triCoo[0])||std::isnan(contactPoint_triCoo[1])||std::isnan(contactPoint_triCoo[2]))
                      {
                          error->warning("Velocity calculation for moving wall failed. Using zero velocity for damping, which may result in over- or underdamping.");
                          for(int mm=0;mm<3;mm++) vwall[mm]=0.;
                      }
                      else
                      {
                          //to reduce numerical error, normalize the coordinates manually
                          normalize_bary(contactPoint_triCoo);

                          vnode=FixMeshGranList[iFMG]->STLdata->v_node[iTri];
                          for (int mm=0;mm<3;mm++) {
                             vwall[mm]=contactPoint_triCoo[0]*vnode[0][mm]+contactPoint_triCoo[1]*vnode[1][mm]+contactPoint_triCoo[2]*vnode[2][mm];
                          }
                      }
                      
                  }
                  else if(FixMeshGranList[iFMG]->STLdata->conveyor) for(int mm=0;mm<3;mm++) vwall[mm]=FixMeshGranList[iFMG]->STLdata->v_node[iTri][0][mm]; //all nodes have same vel
                  else for(int mm=0;mm<3;mm++) vwall[mm]=0.;

                  atom_type_wall=FixMeshGranList[iFMG]->atom_type_wall;

                  //get the parameters needed to resolve the contact
                  deriveContactModelParams(i,-deltan,kn,kt,gamman,gammat,xmu); //deltan>0 in this function

                  delr=radius[i]+deltan; //deltan<0
                  rsq=delr*delr;
                  dx = -delr*en0[0];
                  dy = -delr*en0[1];
                  dz = -delr*en0[2];

                  shr=NULL;
                  if (shear!=NULL) shr=shear[i][iContactList];

                  if (FixMeshGranList[iFMG]->analyseStress)
                  {
                        //calculate force on particle and wall triangle
                        vectorCopy3D(f[i],wallforce);
                        compute_force(i,rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],radius[i],rmass[i],shr,kn,kt,gamman,gammat,xmu);
                        vectorSubtract3D(wallforce,f[i],wallforce);
                        FixMeshGranList[iFMG]->add_particle_contribution(wallforce,contactPoint,iTri);
                  }
                  else compute_force(i,rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],radius[i],rmass[i],shr,kn,kt,gamman,gammat,xmu);

              }
              
              else if(npartners[i]) delflag[j] = 1;
          }

          for(int j=0; j < nTriList[i]; j++)
          {
              //unset touching neighs
              if(delflag[j]) remove_from_contact_list(i,tri_neighlist[i][j][0]/*iFMG*/,tri_neighlist[i][j][1]/*iTri*/);
          }
      }
    }
  }
  
  int maxpartners_all;
  MPI_Allreduce(&maxpartners,&maxpartners_all,1,MPI_INT,MPI_MAX,world);
  
  while(maxpartners<maxpartners_all)
  {
      grow_arrays_maxtritouch(atom->nmax);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   add tri to contact list; if already there, just return the index
   if contact list not large enough, grow it
------------------------------------------------------------------------- */

inline int FixWallGranHookeHistory::add_to_contact_list(int i,int iFMG,int iTri)
{
    for(int k=0;k<npartners[i];k++)
        if(partner[i][k][0]==iFMG && partner[i][k][1]==iTri) return k;

    if(npartners[i]==maxpartners) grow_arrays_maxtritouch(atom->nmax);

    partner[i][npartners[i]][0]=iFMG;
    partner[i][npartners[i]][1]=iTri;
    npartners[i]++;

    shear_transition(i,iFMG,iTri);

    return (npartners[i]-1);
}

inline void FixWallGranHookeHistory::shear_transition(int i,int iFMG,int iTri)
{
    
    for(int j = 0; j < (npartners[i]-1); j++)
    {
        if (partner[i][j][0] != iFMG) continue;
        int iPartner = partner[i][j][1];
        if(FixMeshGranList[iFMG]->STLdata->are_coplanar_neighs(iTri,iPartner))
        {
            
            double shrmagsqr = shear[i][j][0]*shear[i][j][0]+shear[i][j][1]*shear[i][j][1]+shear[i][j][2]*shear[i][j][2];
            if(shrmagsqr > 0.)
            {
                shear[i][npartners[i]-1][0] = shear[i][j][0];
                shear[i][npartners[i]-1][1] = shear[i][j][1];
                shear[i][npartners[i]-1][2] = shear[i][j][2];
                
            }
        }
    }
}

/* ----------------------------------------------------------------------
   remove tri from contact list; also reset shear history
------------------------------------------------------------------------- */

inline void FixWallGranHookeHistory::remove_from_contact_list(int i,int iFMG,int iTri)
{
    for(int k=0;k<npartners[i];k++)
    {
        if(partner[i][k][0]==iFMG && partner[i][k][1]==iTri)
        {
            if(npartners[i]>1) //swap with last contact if more than 1
            {
                partner[i][k][0] = partner[i][npartners[i]-1][0];
                partner[i][k][1] = partner[i][npartners[i]-1][1];
                shear[i][k][0] = shear[i][npartners[i]-1][0];
                shear[i][k][1] = shear[i][npartners[i]-1][1];
                shear[i][k][2] = shear[i][npartners[i]-1][2];
            }

            partner[i][npartners[i]-1][0] = -1;
            partner[i][npartners[i]-1][1] = -1;
            resetShearHistory(i,npartners[i]-1);
            npartners[i]--;
        }
    }
    return;
}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHookeHistory::resetShearHistory(int i,int j)
{
    if (this->shear!=NULL){
        shear[i][j][0] = 0.0;
	    shear[i][j][1] = 0.0;
	    shear[i][j][2] = 0.0;
    }
    else error->all("wall/gran: shear history error");
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::compute_force(int ip, double rsq, double dx, double dy, double dz,
				double *vwall, double *v,
				double *f, double *omega, double *torque,
				double radius, double mass, double *shear,
				double kn,double kt,double gamman, double gammat, double xmu)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,rinv,rsqinv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = radius*omega[0] * rinv;
  wr2 = radius*omega[1] * rinv;
  wr3 = radius*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = gamman*vnnr*rsqinv;       
  ccel = kn*(radius-r)*rinv - damp;
  
  if(cohesionflag)
  {
      double Fn_coh;
      addCohesionForce(ip, r, Fn_coh);
      ccel-=Fn_coh*rinv;
  }

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  shear[0] += vtr1*dt;
  shear[1] += vtr2*dt;
  shear[2] += vtr3*dt;
  shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] + shear[2]*shear[2]);

  // rotate shear displacements

  rsht = shear[0]*dx + shear[1]*dy + shear[2]*dz;
  rsht = rsht*rsqinv;
  shear[0] -= rsht*dx;
  shear[1] -= rsht*dy;
  shear[2] -= rsht*dz;

  // tangential forces = shear + tangential velocity damping

  fs1 = - (kt*shear[0] + gammat*vtr1);  
  fs2 = - (kt*shear[1] + gammat*vtr2);  
  fs3 = - (kt*shear[2] + gammat*vtr3);  

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  if (fs > fn) {
    if (shrmag != 0.0) {
      shear[0] = (fn/fs) * (shear[0] + gammat*vtr1/kt) -      
	gammat*vtr1/kt;
      shear[1] = (fn/fs) * (shear[1] + gammat*vtr2/kt) -      
	gammat*vtr2/kt;
      shear[2] = (fn/fs) * (shear[2] + gammat*vtr3/kt) -      
	gammat*vtr3/kt;
      fs1 *= fn/fs ;
      fs2 *= fn/fs;
      fs3 *= fn/fs;
    } else fs1 = fs2 = fs3 = 0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  torque[0] -= radius*tor1;
  torque[1] -= radius*tor2;
  torque[2] -= radius*tor3;
}

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void FixWallGranHookeHistory::addCohesionForce(int &ip, double &r, double &Fn_coh) 
{
    //r is the distance between the sphere center and wall
    double Acont = (reff_wall*reff_wall-r*r)*M_PI; //contact area sphere-wall
    Fn_coh=mpg->cohEnergyDens[itype][atom_type_wall]*Acont;
}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHookeHistory::deriveContactModelParams(int &ip, double deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu) 
{
    
    double meff_wall = atom->rmass[ip];
    if(fr&&fr->body[ip]>=0)
    {
        meff_wall=fr->masstotal[fr->body[ip]];  
    }

    kn=16./15.*sqrt(reff_wall)*mpg->Yeff[itype][atom_type_wall]*pow(15.*meff_wall*pow(mpg->charVel,2.)/(16.*sqrt(reff_wall)*mpg->Yeff[itype][atom_type_wall]),0.2);
    kt=kn;
    gamman=sqrt(4.*meff_wall*kn/(1.+(M_PI/mpg->coeffRestLog[itype][atom_type_wall])*(M_PI/mpg->coeffRestLog[itype][atom_type_wall])));
    gammat=gamman;
    xmu=mpg->coeffFrict[itype][atom_type_wall];

    if (dampflag == 0) gammat = 0.0;

    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

/* ----------------------------------------------------------------------
   zero wall forces
------------------------------------------------------------------------- */

inline void FixWallGranHookeHistory::reset_wall_forces()
{
      for (int i=0;i<nFixMeshGran;i++)
      {
          for(int j=0;j<FixMeshGranList[i]->nTri;j++)
          {
              for(int k=0;k<3;k++) FixMeshGranList[i]->STLdata->f_tri[j][k]=0.;
          }
      }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixWallGranHookeHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += maxpartners * 2*nmax * sizeof(int);

  if(!shearhistory)return bytes;
  bytes += maxpartners * 3*nmax * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::grow_arrays(int nmax)
{
  npartners = (int*)(memory->srealloc(npartners, nmax*sizeof(int), "fix_wall_gran:npartners"));
  partner = memory->grow_3d_int_array(partner,nmax,maxpartners,2,"fix_wall_gran:partner");

  if(!shearhistory)return;
  shear   = memory->grow_3d_double_array(shear,nmax,maxpartners,3,"fix_wall_gran:shear");

}

/* ----------------------------------------------------------------------
   grow max. # of touching tris 
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::grow_arrays_maxtritouch(int nmax)
{
  if(comm->me==0) {
    if(screen) fprintf(screen,"INFO: Maxmimum number of particle-tri contacts >%d at step %d, growing array...",maxpartners,update->ntimestep);
    if(logfile) fprintf(logfile,"INFO: Maxmimum number of particle-tri contacts >%d at step %d, growing array...",maxpartners,update->ntimestep);
  }
  maxpartners+=DELTA_TRI_CONTACTS;
  
  int ***partner_g = memory->create_3d_int_array(nmax,maxpartners,2,"fix_wall_gran:partner_g");

  for(int i=0;i<nmax;i++)
    for(int j=0;j<maxpartners-DELTA_TRI_CONTACTS;j++)
      for(int k=0;k<2;k++)
         partner_g[i][j][k]=partner[i][j][k];

  int ***h = partner;
  partner=partner_g;
  memory->destroy_3d_int_array(h);

  if(shearhistory)
  {
      double ***shear_g = memory->create_3d_double_array(nmax,maxpartners,3,"fix_wall_gran:shear_g");

      for(int i=0;i<nmax;i++)
        for(int j=0;j<maxpartners-DELTA_TRI_CONTACTS;j++)
          for(int k=0;k<3;k++)
             shear_g[i][j][k]=shear[i][j][k];

      double ***h2=shear;
      shear=shear_g;
      memory->destroy_3d_double_array(h2);
  }

  if(comm->me==0) {
      if(screen) fprintf(screen,"done!\n");
      if(logfile) fprintf(logfile,"done!\n");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::copy_arrays(int i, int j)
{
  npartners[j]=npartners[i];
  for (int k=0;k<maxpartners;k++)
  {
     partner[j][k][0] = partner[i][k][0];
     partner[j][k][1] = partner[i][k][1];
  }

  if(!shearhistory)return;

  for (int k=0;k<maxpartners;k++){
    shear[j][k][0] = shear[i][k][0];
    shear[j][k][1] = shear[i][k][1];
    shear[j][k][2] = shear[i][k][2];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::set_arrays(int i)
{
  npartners[i]=0;
  for (int k=0;k<maxpartners;k++)
  {
     partner[i][k][0] = -1;
     partner[i][k][1] = -1;
  }

  if(!shearhistory) return;
  for (int k=0;k<maxpartners;k++) shear[i][k][0] = shear[i][k][1] = shear[i][k][2] = 0.0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixWallGranHookeHistory::pack_exchange(int i, double *buf)
{
  int m=0;

  buf[m++]=static_cast<int>(npartners[i]);
  for (int k=0;k<maxpartners;k++)
  {
     buf[m++] = static_cast<int>(partner[i][k][0]);
     buf[m++] = static_cast<int>(partner[i][k][1]);
  }

  if(!shearhistory) return m;

  for (int k=0;k<maxpartners;k++){
    buf[m++] = shear[i][k][0];
    buf[m++] = shear[i][k][1];
    buf[m++] = shear[i][k][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixWallGranHookeHistory::unpack_exchange(int nlocal, double *buf)
{
  int m=0;

  npartners[nlocal]=static_cast<int>(buf[m++]);

  for (int k=0;k<maxpartners;k++)
  {
     partner[nlocal][k][0] = static_cast<int>(buf[m++]);
     partner[nlocal][k][1] = static_cast<int>(buf[m++]);
  }

  if(!shearhistory)return m;

  for (int k=0;k<maxpartners;k++){
    shear[nlocal][k][0] = buf[m++];
    shear[nlocal][k][1] = buf[m++];
    shear[nlocal][k][2] = buf[m++];
  }
  return m;
}

/* ----------------------------------------------------------------------
   register and restart fix TriNeighlist now,               
   because otherwise it would be too late to restart it
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::write_restart(FILE * fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<double>(maxpartners);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}
/* ----------------------------------------------------------------------*/
void FixWallGranHookeHistory::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  int maxpartners_restart = static_cast<int> (list[n++]);

  if(maxpartners_restart<1) error->all("Error re-starting fix wall/gran, corrupt restart file?");

  while(maxpartners<maxpartners_restart)
  {
      if(comm->me==0) fprintf(screen,"INFO: Growing number of particle-tri contacts out of restart\n");
      grow_arrays_maxtritouch(atom->nmax);
  }

  if ((wallstyle==MESHGRAN_FWGHH) && (fix_tri_neighlist==NULL)) registerTriNeighlist(0);
  if (wallstyle==MESHGRAN_FWGHH) fix_tri_neighlist->do_warn_dangerous=0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixWallGranHookeHistory::pack_restart(int i, double *buf)
{
  int m = 0;

  buf[m++] = 1+1+maxpartners*2;

  buf[m++]=static_cast<int>(npartners[i]);
  for (int k=0;k<maxpartners;k++)
  {
     buf[m++] = static_cast<int>(partner[i][k][0]);
     buf[m++] = static_cast<int>(partner[i][k][1]);
  }

  if(!shearhistory) return 0;

  buf[0] += maxpartners*3;
  for (int k=0;k<maxpartners;k++){
    buf[m++] = shear[i][k][0];
    buf[m++] = shear[i][k][1];
    buf[m++] = shear[i][k][2];
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixWallGranHookeHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartners[nlocal] = static_cast<int>(extra[nlocal][m++]);

  if(npartners[nlocal]>maxpartners)  error->all("Internal error: Failed to grow maxpartners, incompatible restart file.");

  for (int k=0;k<maxpartners;k++){
    partner[nlocal][k][0] = static_cast<int>(extra[nlocal][m++]);
    partner[nlocal][k][1] = static_cast<int>(extra[nlocal][m++]);
  }

  if(!shearhistory) return;

  for (int k=0;k<maxpartners;k++){
    shear[nlocal][k][0] = extra[nlocal][m++];
    shear[nlocal][k][1] = extra[nlocal][m++];
    shear[nlocal][k][2] = extra[nlocal][m++];
  }

}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixWallGranHookeHistory::maxsize_restart()
{
  if(!shearhistory) return 1+1+maxpartners*2;
  return 1+1+maxpartners*5;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixWallGranHookeHistory::size_restart(int nlocal)
{
  if(!shearhistory) return 1+1+maxpartners*2;
  return 1+1+maxpartners*5;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistory::reset_dt()
{
  dt = update->dt;
}
