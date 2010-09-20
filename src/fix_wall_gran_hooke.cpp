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
#include "fix_wall_gran_hooke.h"
#include "pair_gran_hooke.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"

using namespace LAMMPS_NS;

enum{MESHGRAN,XPLANE,YPLANE,ZPLANE,ZCYLINDER};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixWallGranHooke::FixWallGranHooke(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistory(lmp, narg, arg)
{
    shearhistory=0;
    initSubstyle();
}

/* ---------------------------------------------------------------------- */

int FixWallGranHooke::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHooke::initSubstyle()
{
  if (!(force->pair_match("gran/hooke",0))) error->all("Fix wall/gran/hooke can only be used together with pair style gran/hooke");

  // initialize as if particle is not touching wall
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
  {
      npartners[i]=0;
      for (int k=0;k<maxpartners;k++)
      {
          partner[i][k][0] = partner[i][k][0] = -1;
      }
  }

}

/* ---------------------------------------------------------------------- */

void FixWallGranHooke::init()
{
  dt = update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if ((wallstyle==MESHGRAN) && (fix_tri_neighlist==NULL)) registerTriNeighlist(1);

  if (!(force->pair_match("gran/hooke",0))) error->all("Fix wall/gran/hooke can only be used together with pair style gran/hooke");

  mpg = static_cast<PairGranHooke*>(force->pair)->getMatProp();

  fr=NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
  {
      if(strncmp(modify->fix[ifix]->style,"rigid",5)==0) fr=static_cast<FixRigid*>(modify->fix[ifix]);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranHooke::resetShearHistory(int i,int j)
{
}

/* ---------------------------------------------------------------------- */

inline int FixWallGranHooke::add_to_contact_list(int i,int iFMG,int iTri)
{
    for(int k=0;k<npartners[i];k++)
        if(partner[i][k][0]==iFMG && partner[i][k][1]==iTri) return k;

    if(npartners[i]==maxpartners) grow_arrays_maxtritouch(atom->nmax);

    partner[i][npartners[i]][0]=iFMG;
    partner[i][npartners[i]][1]=iTri;
    npartners[i]++;

    //skip shear transition here

    return (npartners[i]-1);
}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHooke::remove_from_contact_list(int i,int iFMG,int iTri)
{
    for(int k=0;k<npartners[i];k++)
        if(partner[i][k][0]==iFMG && partner[i][k][0]==iTri)
        {
            partner[i][k][0] = partner[i][npartners[i]-1][0];
            partner[i][k][1] = partner[i][npartners[i]-1][1];

            partner[i][npartners[i]-1][0] = -1;
            partner[i][npartners[i]-1][1] = -1;
            npartners[i]--;
        }
    return;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHooke::compute_force(int ip, double rsq, double dx, double dy, double dz,
			double *vwall, double *v,
			double *f, double *omega, double *torque,
			double radius, double mass, double *shear, /*shear is not used!*/
			double kn,double kt,double gamman, double gammat, double xmu)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3,rinv,rsqinv;

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

  meff = mass;
  if(fr&&fr->body[ip]>=0) meff=fr->masstotal[fr->body[ip]];  
  damp = gamman*vnnr*rsqinv;   
  ccel = kn*(radius-r)*rinv - damp;
  
  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = gammat*vrel;         
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // tangential force due to tangential velocity damping

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

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

