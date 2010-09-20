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
   A poor man's triangle neighbor list
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_tri_neighlist.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "comm.h"
#include "fix_wall_gran_hooke_history.h"
#include "triSpherePenetration.h"
#include "fix_propertyPerAtom.h"
#include "fix_meshGran.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20
#define SMALL_DELTA skin/(70.*M_PI)

#define DELTA_TRI_LIST 10 //number of touching triangles that a particle can have [grow size]

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#include "debug_single_CAD.h"

/* ---------------------------------------------------------------------- */

FixTriNeighlist::FixTriNeighlist(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    if (narg < 4) error->all("Illegal fix neighlist/tri command");

    caller_id=new char[strlen(arg[3]) + 1];
    strcpy(caller_id,arg[3]);

    int ifix=modify->find_fix(caller_id);
    if (ifix<0) error->all("Illegal fix id provided to neighlist/tri command");

    if(strncmp(modify->fix[ifix]->style,"wall/gran",9)!=0) error->all("Illegal fix id provided to neighlist/tri command (is not of type wall/gran)");

    caller=static_cast<FixWallGranHookeHistory*>(modify->fix[ifix]);
    nFixMeshGran=static_cast<FixWallGranHookeHistory*>(caller)->nFixMeshGran;
    FixMeshGranList=static_cast<FixWallGranHookeHistory*>(caller)->FixMeshGranList;

    //initialise members
    bsubboxlo=new double[3];
    bsubboxhi=new double[3];
    buildNeighList=0;
    nTriList=NULL;
    tri_neighlist=NULL;
    delflag=NULL;

    create_attribute = 1;

    // perform initial allocation of atom-based arrays
    // register with Atom class

    maxwalllist=DELTA_TRI_LIST;
    grow_arrays(atom->nmax);
    atom->add_callback(0);
    atom->add_callback(1);

    // initialize with zero

    int nmax = atom->nmax;
    for (int i = 0; i < nmax; i++){
        nTriList[i]=0;
        for (int k=0;k<maxwalllist;k++){
          tri_neighlist[i][k][0] = 0;
          tri_neighlist[i][k][1] = 0;
        }
    }

    force_reneighbor = 1;
    next_reneighbor = update->ntimestep + 1;

    do_warn=1;
    do_warn_dangerous=1; 
}

/* ---------------------------------------------------------------------- */

FixTriNeighlist::~FixTriNeighlist()
{
    delete []nTriList;
    memory->destroy_3d_int_array(tri_neighlist);

    // unregister callbacks to this fix from Atom class
    atom->delete_callback(id,0);
    atom->delete_callback(id,1);
}

/* ---------------------------------------------------------------------- */
void FixTriNeighlist::init()
{
    
    next_reneighbor = update->ntimestep + 1;
    
    do_warn=1;
    buildNeighList=1;
}

/* ---------------------------------------------------------------------- */

int FixTriNeighlist::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTriNeighlist::pre_neighbor()
{
    
    buildNeighList=1;

    FixMeshGran* fmg;
    for(int iList=0;iList<nFixMeshGran;iList++)
    {
        fmg=FixMeshGranList[iList];
        
        if (fmg->STLdata->movingMesh) fmg->STLdata->before_rebuild();
    }

}

/* ---------------------------------------------------------------------- */

void FixTriNeighlist::pre_force(int vflag)
{
    int iList;
    FixMeshGran* fmg;
    double **vatom=atom->v;
    double dt=update->dt;
    double maxmoveSq;
    bool excl;

    if (!buildNeighList && decide_rebuild()) {
        next_reneighbor = update->ntimestep + 1;
        
    }

    //if neighbor list was not built this step, return now
    if (!buildNeighList) return;

    buildNeighList=0;

    //unset non-touching contacts
    unset_nontouching();

    int neigh_style=neighbor->style;
    if(neigh_style!=1) error->all("Please use style 'bin' in the 'neighbor' command together with triangular walls");
    double skin=neighbor->skin;
    double skin_safety=1.;                    
    double cutneighmax=neighbor->cutneighmax; 
    
    int mbinx=neighbor->mbinx;
    int mbiny=neighbor->mbiny;
    int *bins=neighbor->bins;
    int *binhead=neighbor->binhead;

    //some allocations
    int ix,ixmin,ixmax,iy,iymin,iymax,iz,izmin,izmax,binmin,binmax;
    int iTri,nTri,ibin;
    double *tri_xmin=new double[3];
    double *tri_xmax=new double[3];

    double dist,treshold;

    int dangerous_build=0;

    for(int j=0;j<3;j++)
    {
        bsubboxlo[j] = domain->sublo[j] - comm->cutghost[j];
        bsubboxhi[j] = domain->subhi[j] + comm->cutghost[j];
    }

    double en0[3]; //just a dummy

    flag_old_list();

    int iatom_bin;

    for(iList=0;iList<nFixMeshGran;iList++)
    {
        
        fmg=FixMeshGranList[iList];
        nTri=fmg->nTri;

        if (fmg->STLdata->movingMesh) skin_safety = 2.2;
        else skin_safety=1.1;

        treshold=skin*0.5*skin_safety;

        for(iTri=0;iTri<nTri;iTri++)
        {
            
            fmg->STLdata->getTriBoundBox(iTri,tri_xmin,tri_xmax,(cutneighmax+(skin_safety-1.)*skin)+SMALL_DELTA);

            if (!check_box_overlap(bsubboxlo,bsubboxhi,tri_xmin,tri_xmax)) continue;

            for(int j=0;j<3;j++)
            {
                if(tri_xmin[j]<bsubboxlo[j]) tri_xmin[j]=bsubboxlo[j]+SMALL_DELTA;
                if(tri_xmax[j]>bsubboxhi[j]) tri_xmax[j]=bsubboxhi[j]-SMALL_DELTA;
            }

            binmin=neighbor->coord2bin(tri_xmin);
            ixmin = (binmin % (mbiny*mbinx)) % mbinx;
            iymin = static_cast<int>(round(static_cast<double>(((binmin - ixmin) % (mbiny*mbinx)))/static_cast<double>(mbinx)));
            izmin = static_cast<int>(round(static_cast<double>((binmin - ixmin - iymin*mbiny))/static_cast<double>((mbiny*mbinx))));
            binmax= neighbor->coord2bin(tri_xmax);
            ixmax = (binmax % (mbiny*mbinx)) % mbinx;
            iymax = static_cast<int>(round(static_cast<double>(((binmax - ixmin) % (mbiny*mbinx)))/static_cast<double>(mbinx)));
            izmax = static_cast<int>(round(static_cast<double>((binmax - ixmin - iymin*mbiny))/static_cast<double>((mbiny*mbinx))));

            ix=ixmin;
            iy=iymin;
            iz=izmin;
            while(ix<=ixmax)
            {
                iy=iymin;
                iz=izmin;
                while(iy<=iymax)
                {

                    iz=izmin;
                    while(iz<=izmax)
                    {

                        ibin=(iz*mbiny*mbinx + iy*mbinx + ix);
                        
                        iatom_bin=binhead[ibin];
                        while(iatom_bin!=-1 && iatom_bin<atom->nlocal)
                        {

                            dist=TRISPHERE::resolveTriSphereContact(lmp,iatom_bin,iTri,fmg->STLdata,F_SHRINKAGE,en0,true,treshold,excl); //call it with shrinkage =0
                            
                            if(dist<treshold) 
                            {
                                
                                int append = addTriToNeighList(iatom_bin,iTri,iList);
                                dangerous_build = MAX(dangerous_build, check_dangerous(dist,append) );

                                if(dangerous_build&&do_warn_dangerous)
                                {
                                    
                                }

                            }
                            if (bins!=NULL) iatom_bin=bins[iatom_bin];
                            else iatom_bin=-1;
                        }
                        iz++;
                    }
                    iy++;
                }
                ix++;
            }
        }
    }
    
    int maxwalllist_all;
    MPI_Allreduce(&maxwalllist,&maxwalllist_all,1,MPI_INT,MPI_MAX,world);
    
    while(maxwalllist<maxwalllist_all)
    {
        grow_arrays_maxtritouch(atom->nmax);
    }

    clear_old_entries();
    delete []tri_xmin;
    delete []tri_xmax;

    check_dangerous_build(dangerous_build);

}

/* ------------------------------------------------------------------------------------------------------- */

inline void FixTriNeighlist::check_dangerous_build(int dangerous_build)
{
    int dangerous_build_all=0;
    MPI_Allreduce(&dangerous_build,&dangerous_build_all,1,MPI_INT,MPI_MAX,world);
    if(do_warn_dangerous && dangerous_build_all)
    {
        if(screen&&comm->me==0) fprintf(screen,"At step %d\n",update->ntimestep);
        error->warning("Dangerous build in triangle neighbor list.");
        if(logfile) fprintf(logfile,"WARNING: Dangerous build in triangle neighbor list.\n");
        
    }
    
    do_warn_dangerous=1;
}

/* ------------------------------------------------------------------------------------------------------- */
inline void FixTriNeighlist::flag_old_list()
{
    //clear neigh list
    for(int i = 0; i < atom->nlocal; i++){
        if(nTriList[i]>maxwalllist) error->one("Internal error: Inconsistent neigh list");
        for (int j=0;j<nTriList[i];j++){
            tri_neighlist[i][j][0]=-(tri_neighlist[i][j][0]+2);
        }
    }
}

/* ------------------------------------------------------------------------------------------------------- */
inline void FixTriNeighlist::clear_old_entries()
{
    //clear old neigh list entries
    int iMesh, iTri;
    for(int i = 0; i < atom->nlocal; i++){
        for (int j=0;j<nTriList[i];j++){
            iMesh=tri_neighlist[i][j][0];
            iTri=tri_neighlist[i][j][1];

            if(iMesh==-1)error->all("Internal error: Inconsistent neigh list");
            while(iMesh<-1)
            {
                tri_neighlist[i][j][0] = tri_neighlist[i][nTriList[i]-1][0];
                tri_neighlist[i][j][1] = tri_neighlist[i][nTriList[i]-1][1];
                tri_neighlist[i][nTriList[i]-1][0] = -1;
                tri_neighlist[i][nTriList[i]-1][1] = -1;
                nTriList[i]--;

                iMesh=tri_neighlist[i][j][0];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------- */

inline int FixTriNeighlist::decide_rebuild()
{
  if (neighbor->ago >= neighbor->delay && neighbor->ago % neighbor->every == 0)
  {
    if (neighbor->build_once) return 0;
    if (neighbor->dist_check == 0) return 1;
    return check_distance();
  }
  else return 0;
}

/* ------------------------------------------------------------------------------------------------------- */

inline int FixTriNeighlist::check_distance()
{
  double delx,dely,delz,rsq;

  double ***nodes_now,***nodes_lastRe;
  FixMeshGran* fmg;

  int flag = 0;
  
  for(int iList=0;iList<nFixMeshGran;iList++)
  {
      fmg=FixMeshGranList[iList];

      if (fmg->STLdata->movingMesh)
      {
          int nTri=fmg->STLdata->nTri;
          nodes_now=fmg->STLdata->node;
          nodes_lastRe=fmg->STLdata->node_lastRe;

          for(int iTri=0;iTri<nTri;iTri++)
          {
              
              for(int j=0;j<3;j++)
              {
                  delx = nodes_now[iTri][j][0] - nodes_lastRe[iTri][j][0];
                  dely = nodes_now[iTri][j][1] - nodes_lastRe[iTri][j][1];
                  delz = nodes_now[iTri][j][2] - nodes_lastRe[iTri][j][2];
                  rsq = delx*delx + dely*dely + delz*delz;

                  if (rsq > neighbor->triggersq) flag = 1;
              }
          }
      }
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  return flagall;
}

/* ------------------------------------------------------------------------------------------------------- */

inline int FixTriNeighlist::addTriToNeighList(int ip,int iTri, int iMeshGran)
{
    
    int ip_nCont=nTriList[ip];

    int jMeshGran,jTri;
    for(int j=0;j<ip_nCont;j++)
    {
        jMeshGran = tri_neighlist[ip][j][0];
        if (jMeshGran<0) jMeshGran = -jMeshGran - 2;
        jTri      = tri_neighlist[ip][j][1];

        if ((iMeshGran==jMeshGran)&&(iTri==jTri))
        {
            
            if(tri_neighlist[ip][j][0]<0) tri_neighlist[ip][j][0] = jMeshGran;
            return 0;
        }
    }

    //   tri has not been found yet in list, append it
    
    if(ip_nCont==maxwalllist)
    {
        
        grow_arrays_maxtritouch(atom->nmax);
        
    }
    if(ip_nCont==maxwalllist) error->one("Grow of number of particle-triangle contacts failed");

    tri_neighlist[ip][ip_nCont][0]=iMeshGran;
    tri_neighlist[ip][ip_nCont][1]=iTri;
    nTriList[ip]++;

    return 1;
}

/* ---------------------------------------------------------------------- */
inline int FixTriNeighlist::check_dangerous(double dist, int append)
{
    
    if(dist>=0. || append!=1 ) return 0;

    return 1;
}

/* ---------------------------------------------------------------------- */

inline void FixTriNeighlist::unset_nontouching()
{
    //get contact list
    int *npartners = caller->npartners;
    int ***partner = caller->partner;
    int iTri,iFMG;
    double deltan;
    double en0[3];
    bool excl;

    for(int i=0;i<atom->nlocal;i++)
    {
        for(int p=0;p<npartners[i];p++)
        {
            iFMG = partner[i][p][0];
            iTri = partner[i][p][1];
            deltan=TRISPHERE::resolveTriSphereContact(lmp,i,iTri,caller->FixMeshGranList[iFMG]->STLdata,F_SHRINKAGE,en0,false,0.,excl);
            if(deltan>=0.||excl)
            {
                caller->remove_from_contact_list(i,iFMG,iTri);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

bool FixTriNeighlist::check_box_overlap(double *xmin1,double *xmax1,double *xmin2,double *xmax2 )
{
    bool overlap=true;

    for (int i=0;i<3;i++)
      overlap = overlap && ( (xmin1[i]>xmin2[i] && xmin1[i]<xmax2[i]) || (xmin2[i]>xmin1[i] && xmin2[i]<xmax1[i]) );
    return overlap;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixTriNeighlist::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += maxwalllist* 3*nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixTriNeighlist::grow_arrays(int nmax)
{
    nTriList = (int*)(lmp->memory->srealloc(nTriList, nmax*sizeof(int), "fix_tri_neighlist:nTriList"));
    tri_neighlist = memory->grow_3d_int_array(tri_neighlist,nmax,maxwalllist,2,"fix_tri_neighlist:tri_neighlist");

    delflag = (int*)(memory->srealloc(delflag, maxwalllist*sizeof(int), "fix_tri_neighlist:delflag"));
}

/* ----------------------------------------------------------------------
   grow 2nd dim of tri neighlist
------------------------------------------------------------------------- */

void FixTriNeighlist::grow_arrays_maxtritouch(int nmax)
{
    if(comm->me==0) {
            if(screen) fprintf(screen,"INFO: Maxmimum number of particle-tri neighbors >%d at step %d, growing array...",maxwalllist,update->ntimestep);
            if(logfile) fprintf(logfile,"INFO: Maxmimum number of particle-tri neighbors >%d at step %d, growing array...",maxwalllist,update->ntimestep);
    }

    maxwalllist+=DELTA_TRI_LIST;
    
    //grow tri_neighlist
    int ***tri_neighlist_g = memory->create_3d_int_array(nmax,maxwalllist,2,"fix_tri_neighlist:tri_neighlist_g");

    for(int i=0;i<nmax;i++)
        for(int j=0;j<maxwalllist-DELTA_TRI_LIST;j++)
            for(int k=0;k<2;k++)
                tri_neighlist_g[i][j][k]=tri_neighlist[i][j][k];

    //clean up
    int ***h=tri_neighlist;
    tri_neighlist=tri_neighlist_g;
    memory->destroy_3d_int_array(h);

    delflag = (int*)(memory->srealloc(delflag, maxwalllist*sizeof(int), "fix_tri_neighlist:delflag"));

    if(comm->me==0) {
        if(screen) fprintf(screen,"done!\n");
        if(logfile) fprintf(logfile,"done!\n");
    }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixTriNeighlist::copy_arrays(int i, int j)
{
  nTriList[j]=nTriList[i];
  for (int k=0;k<maxwalllist;k++){
    tri_neighlist[j][k][0] = tri_neighlist[i][k][0];
    tri_neighlist[j][k][1] = tri_neighlist[i][k][1];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixTriNeighlist::set_arrays(int i)
{
  nTriList[i]=0;
  for (int k=0;k<maxwalllist;k++) tri_neighlist[i][k][0] = tri_neighlist[i][k][1] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixTriNeighlist::pack_exchange(int i, double *buf)
{
  int m=0;
  buf[m++] = static_cast<double>(nTriList[i]);
  for (int k=0;k<maxwalllist;k++){
    buf[m++] = static_cast<double>(tri_neighlist[i][k][0]);
    buf[m++] = static_cast<double>(tri_neighlist[i][k][1]);
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values into local atom-based arrays after exchange
------------------------------------------------------------------------- */

int FixTriNeighlist::unpack_exchange(int nlocal, double *buf)
{
  int m=0;
  nTriList[nlocal]=static_cast<int>(buf[m++]);
  for (int k=0;k<maxwalllist;k++){
    tri_neighlist[nlocal][k][0] = static_cast<int>(buf[m++]);
    tri_neighlist[nlocal][k][1] = static_cast<int>(buf[m++]);
  }
  return m;
}

