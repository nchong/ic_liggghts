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
#include "fix_wall_gran_hertz_history.h"
#include "pair_gran_hertz_history.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "mech_param_gran.h"
#include "fix_rigid.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

enum{MESHGRAN,XPLANE,YPLANE,ZPLANE,ZCYLINDER};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

/* ---------------------------------------------------------------------- */

FixWallGranHertzHistory::FixWallGranHertzHistory(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistory(lmp, narg, arg)
{
  initSubstyle();
}

/* ---------------------------------------------------------------------- */

int FixWallGranHertzHistory::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHertzHistory::init()
{
  dt = update->dt;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if ((wallstyle==MESHGRAN) && (fix_tri_neighlist==NULL)) registerTriNeighlist(1);

  if (!force->pair_match("gran/hertz/history",0)) error->all("Fix wall/gran/hertz/history can only be used together with pair style gran/hertz/history");

  mpg = static_cast<PairGranHertzHistory*>(force->pair)->getMatProp();

  fr=NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
  {
      if(strncmp(modify->fix[ifix]->style,"rigid",5)==0) fr=static_cast<FixRigid*>(modify->fix[ifix]);
  }
}

void FixWallGranHertzHistory::initSubstyle()
{
  if(!(force->pair_match("gran/hertz/history",0))) error->all("Fix wall/gran/hertz/history can only be used together with pair style gran/hertz/history");

  mpg = static_cast<PairGranHertzHistory*>(force->pair)->getMatProp();
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void FixWallGranHertzHistory::deriveContactModelParams(int &ip, double deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu)  
{

    double meff_wall = atom->rmass[ip];
    if(fr&&fr->body[ip]>=0)
    {
        meff_wall=fr->masstotal[fr->body[ip]];  
    }

    double sqrtval = sqrt(reff_wall*deltan);

    double Sn=2.*(mpg->Yeff[itype][atom_type_wall])*sqrtval;
    double St=8.*(mpg->Geff[itype][atom_type_wall])*sqrtval;

    kn=4./3.*mpg->Yeff[itype][atom_type_wall]*sqrtval;
    kt=St;

    gamman=-2.*sqrtFiveOverSix*mpg->betaeff[itype][atom_type_wall]*sqrt(Sn*meff_wall);
    gammat=-2.*sqrtFiveOverSix*mpg->betaeff[itype][atom_type_wall]*sqrt(St*meff_wall);

    xmu=mpg->coeffFrict[itype][atom_type_wall];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

