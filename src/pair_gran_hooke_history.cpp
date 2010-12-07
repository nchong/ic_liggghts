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

//#define EMIT_DIAGNOSTICS //< Write to [diagnostics.txt]
//#define EMIT_PAIRWISE //< Write to [pairwise_data.csv] for [play/hertz.cu]

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_history.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_pour.h"
#include "fix_shear_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "mech_param_gran.h"
#include "fix_rigid.h"
#include "fix_pour.h"
#include "fix_pour_dev.h"
#include "fix_particledistribution_discrete.h"
#include "fix_pour_legacy.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::PairGranHookeHistory(LAMMPS *lmp) : Pair(lmp)
{
  
  if (!(force->pair_match("gran/hooke", 0)) && !(force->pair_match("gran/hertz", 0))) return;

  single_enable = 0;
  no_virial_compute = 1;
  history = 1;
  fix_history = NULL;

  mpg=new MechParamGran(lmp);
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::~PairGranHookeHistory()
{
  if (fix_history) modify->delete_fix("SHEAR_HISTORY");

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
  }
  delete mpg;

  delete [] onerad_dynamic;
  delete [] onerad_frozen;
  delete [] maxrad_dynamic;
  delete [] maxrad_frozen;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void PairGranHookeHistory::addCohesionForce(int &ip, int &jp,double &r, double &Fn_coh) 
{
    //r is the distance between the sphere's centeres
    double Acont = - M_PI/4 * ( (r-ri-rj)*(r+ri-rj)*(r-ri+rj)*(r+ri+rj) )/(r*r); //contact area of the two spheres
    
    Fn_coh=mpg->cohEnergyDens[itype][jtype]*Acont;
}

/* ---------------------------------------------------------------------- */

inline void PairGranHookeHistory::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu) 
{
    double reff=ri*rj/(ri+rj);
    kn=16./15.*sqrt(reff)*(mpg->Yeff[itype][jtype])*pow(15.*meff*mpg->charVel*mpg->charVel/(16.*sqrt(reff)*mpg->Yeff[itype][jtype]),0.2);
    kt=kn;
    gamman=sqrt(4.*meff*kn/(1.+(M_PI/mpg->coeffRestLog[itype][jtype])*(M_PI/mpg->coeffRestLog[itype][jtype])));
    gammat=gamman;
    xmu=mpg->coeffFrict[itype][jtype];

    if (dampflag == 0) gammat = 0.0;

    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::compute(int eflag, int vflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3;
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  // loop over neighbors of my atoms

#ifdef EMIT_PAIRWISE
  static bool first_call = true;
  FILE *ofile = fopen("pairwise_data.csv", "a");
#endif
#ifdef EMIT_DIAGNOSTICS
  FILE *hist_file = fopen("numneigh.txt", "a");
  int num_contacts_hist[64];
  memset(num_contacts_hist, 0, 64 * sizeof(int));

  static int prev_num_contacts = -1;
  int num_contacts = 0;
  static int step = -1;
  step++;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jnum = numneigh[i];
    num_contacts += jnum;
    num_contacts_hist[jnum]++;
  }
  fprintf(hist_file, "-----\n");
  for (int i=0; i<64; i++) {
    fprintf(hist_file, "%d, %d\n", i, num_contacts_hist[i]);
  }
  fflush(hist_file);
  fclose(hist_file);

  if (num_contacts != prev_num_contacts) {
    prev_num_contacts = num_contacts;
    FILE *num_contacts_file = fopen("numcontacts.txt", "a");
    fprintf(num_contacts_file, "%d, %d\n", step, num_contacts);
    fflush(num_contacts_file);
    fclose(num_contacts_file);
  }
#endif

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

	// unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[3*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
#ifdef EMIT_PAIRWISE
        if (first_call) {
          first_call = false;
          printf("dampflag=%d\n", dampflag);
          printf("eflag=%d vflag=%d cohesionflag=%d\n", eflag, vflag, cohesionflag);
          printf("freeze_group_bit[i]=%d freeze_group_bit[j]=%d\n",
            mask[i] & freeze_group_bit,
            mask[j] & freeze_group_bit);

          // fixed parameters need to be copy+pasted
          printf("float dt = %.16f;\n", dt);
          printf("float Yeff[%d][%d];\n",
            (mpg->max_type+1),
            (mpg->max_type+1));
          printf("float Geff[%d][%d];\n",
            (mpg->max_type+1),
            (mpg->max_type+1));
          printf("float betaeff[%d][%d];\n",
            (mpg->max_type+1),
            (mpg->max_type+1));
          printf("float coeffFrict[%d][%d];\n",
            (mpg->max_type+1),
            (mpg->max_type+1));
          // nb: type 0 is not defined
          for (int p=1; p<(mpg->max_type+1); p++) {
            for (int q=1; q<(mpg->max_type+1); q++) {
              printf("Yeff[%d][%d] = %.16f;\n", p,q,mpg->Yeff[p][q]);
              printf("Geff[%d][%d] = %.16f;\n", p,q,mpg->Geff[p][q]);
              printf("betaeff[%d][%d] = %.16f;\n", p,q,mpg->betaeff[p][q]);
              printf("coeffFrict[%d][%d] = %.16f;\n", p,q,mpg->coeffFrict[p][q]);
            }
          }
          printf("float nktv2p = %.16f;\n", force->nktv2p);
          printf("int typei = %d;\nint typej = %d;\n", type[i], type[j]);
          printf("int max_type = %d;\n", (mpg->max_type+1));
        }
        fprintf(ofile, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f, ",
          x[i][0], x[i][1], x[i][2],
          x[j][0], x[j][1], x[j][2]
        );
        fprintf(ofile, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f, ",
          v[i][0], v[i][1], v[i][2],
          v[j][0], v[j][1], v[j][2]
        );
        fprintf(ofile, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f, ",
          omega[i][0], omega[i][1], omega[i][2],
          omega[j][0], omega[j][1], omega[j][2]
        );
        fprintf(ofile, "%.16f, %.16f, ",
          radius[i], radius[j]);
        fprintf(ofile, "%.16f, %.16f, ",
          rmass[i], rmass[j]);
        //fprintf(ofile, "%.16f, %.16f, ",
        //  mass[type[i]], mass[type[j]]);

        double *shear = &allshear[3*jj];
        fprintf(ofile, "%.16f, %.16f, %.16f, ",
          shear[0], shear[1], shear[2]);
        fprintf(ofile, "%.16f, %.16f, %.16f, ",
          torque[i][0], torque[i][1], torque[i][2]);
        fprintf(ofile, "%.16f, %.16f, %.16f, ",
          f[i][0], f[i][1], f[i][2]);
#endif

	r = sqrt(rsq);
	rinv = 1.0/r;
	rsqinv = 1.0/rsq;

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	// normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr * rsqinv;
	vn2 = dely*vnnr * rsqinv;
	vn3 = delz*vnnr * rsqinv;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	// relative rotational velocity

	wr1 = (radi*omega[i][0] + radj*omega[j][0]) * rinv;
	wr2 = (radi*omega[i][1] + radj*omega[j][1]) * rinv;
	wr3 = (radi*omega[i][2] + radj*omega[j][2]) * rinv;

	// normal forces = Hookian contact + normal velocity damping

	double mi,mj;
	if (rmass) {
	  mi=rmass[i];
	  mj=rmass[j];
	} else {
	  itype = type[i];
	  jtype = type[j];
	  mi=mass[itype];
	  mj=mass[jtype];
	}
	if (fr) //deadcode?
	{
	   if(fr->body[i]>=0) double mi=fr->masstotal[fr->body[i]];  
	   if(fr->body[j]>=0) double mj=fr->masstotal[fr->body[j]];  
	}
	meff=mi*mj/(mi+mj);
	if (mask[i] & freeze_group_bit) meff = mj;
	if (mask[j] & freeze_group_bit) meff = mi;

    double deltan=radsum-r;
	deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu);	 //modified C.K

	damp = gamman*vnnr*rsqinv;  
	ccel = kn*(radsum-r)*rinv - damp;
	
	if (cohesionflag) { 
	    addCohesionForce(i,j,r,Fn_coh);
	    ccel-=Fn_coh*rinv;
	}

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3; //deadcode?
	vrel = sqrt(vrel);

	// shear history effects

	touch[jj] = 1;
	shear = &allshear[3*jj];
        shear[0] += vtr1*dt;
        shear[1] += vtr2*dt;
        shear[2] += vtr3*dt;
        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +
		      shear[2]*shear[2]);

	// rotate shear displacements

	rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
	rsht *= rsqinv;
        shear[0] -= rsht*delx;
        shear[1] -= rsht*dely;
        shear[2] -= rsht*delz;

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
	    fs1 *= fn/fs;
	    fs2 *= fn/fs;
	    fs3 *= fn/fs;
	  } else fs1 = fs2 = fs3 = 0.0;
	}

	// forces & torques

	fx = delx*ccel + fs1;
	fy = dely*ccel + fs2;
	fz = delz*ccel + fs3;
	f[i][0] += fx;
	f[i][1] += fy;
	f[i][2] += fz;

	tor1 = rinv * (dely*fs3 - delz*fs2);
	tor2 = rinv * (delz*fs1 - delx*fs3);
	tor3 = rinv * (delx*fs2 - dely*fs1);
	torque[i][0] -= radi*tor1;
	torque[i][1] -= radi*tor2;
	torque[i][2] -= radi*tor3;

#ifdef EMIT_PAIRWISE
  fprintf(ofile, "%.16f, %.16f, %.16f, ",
    shear[0], shear[1], shear[2]);
  fprintf(ofile, "%.16f, %.16f, %.16f, ",
    torque[i][0], torque[i][1], torque[i][2]);
  fprintf(ofile, "%.16f, %.16f, %.16f\n",
    f[i][0], f[i][1], f[i][2]);
#endif

	if (j < nlocal) {
	  f[j][0] -= fx;
	  f[j][1] -= fy;
	  f[j][2] -= fz;
	  torque[j][0] -= radj*tor1;
	  torque[j][1] -= radj*tor2;
	  torque[j][2] -= radj*tor3;
	}

	if (evflag) ev_tally_xyz(i,j,nlocal,0,
				 0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

#ifdef EMIT_PAIRWISE
  fflush(ofile);
  fclose(ofile);
#endif
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGranHookeHistory::allocate()
{
  allocated = 1;
  int n = atom->ntypes; //ensured elsewhere that this is high enough

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeHistory::settings(int narg, char **arg) 
{
  if (narg != 2) error->all("Illegal pair_style command");

  dampflag = force->inumeric(arg[0]);
  cohesionflag = force->inumeric(arg[1]);

  if (dampflag < 0 || dampflag > 1 || cohesionflag < 0 || cohesionflag > 1)
    error->all("Illegal pair_style command");

  if(cohesionflag && domain->dimension!=3) error->all("Cohesion model valid for 3d simulations only");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGranHookeHistory::coeff(int narg, char **arg)
{
  if (narg > 2) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistory::init_substyle()
{
  mpg->getMaterialParams(1,cohesionflag);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGranHookeHistory::init_style()
{
  int i;

  // error and warning checks

  if(strcmp(update->unit_style,"metal")==0 || strcmp(update->unit_style,"real")==0) error->all("Can not use a non-consistent unit system with pair gran. Please use si,cgs or lj.");

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Pair granular requires atom attributes radius, omega, torque");
  if (comm->ghost_velocity == 0)
    error->all("Pair granular requires ghost atoms store velocity");

  fr=NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
  {
      if(strncmp(modify->fix[ifix]->style,"rigid",5)==0) fr=static_cast<FixRigid*>(modify->fix[ifix]);
  }

  init_substyle();
  
  // need a half neigh list and optionally a granular history neigh list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  if (history) {
    irequest = neighbor->request(this);
    neighbor->requests[irequest]->id = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->granhistory = 1;
    neighbor->requests[irequest]->dnum = 3;
  }

  dt = update->dt;

  // if shear history is stored:
  // check if newton flag is valid
  // if first init, create Fix needed for storing shear history

  if (history && force->newton_pair == 1)
    error->all("Pair granular with shear history requires newton pair off");

  if (history && fix_history == NULL) {
    char **fixarg = new char*[3];
    fixarg[0] = (char *) "SHEAR_HISTORY";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "SHEAR_HISTORY";
    modify->add_fix(3,fixarg);
    delete [] fixarg;
    fix_history = (FixShearHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
  }

  // check for freeze Fix and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // set maxrad_dynamic and maxrad_frozen for each type
  for (i = 1; i <= atom->ntypes; i++)
  onerad_dynamic[i] = onerad_frozen[i] = 0.0;

  // include future Fix pour particles as dynamic

  for (i = 0; i < modify->nfix; i++){
    for(int j=1;j<=atom->ntypes;j++)
    {
        int pour_type = 0;
        double pour_maxrad = 0.0;
        pour_type = j;
        pour_maxrad = modify->fix[i]->max_rad(pour_type);
        onerad_dynamic[pour_type] = MAX(onerad_dynamic[pour_type],pour_maxrad);
    }
  }

  //further dynamic and frozen

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
		MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
		MPI_DOUBLE,MPI_MAX,world);

}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   optional granular history list
------------------------------------------------------------------------- */

void PairGranHookeHistory::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listgranhistory = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGranHookeHistory::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
  cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
  cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);
  return cutoff;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGranHookeHistory::write_restart(FILE *fp)
{
  write_restart_settings(fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGranHookeHistory::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairGranHookeHistory::write_restart_settings(FILE *fp) 
{
  fwrite(&dampflag,sizeof(int),1,fp);
  fwrite(&cohesionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts 
------------------------------------------------------------------------- */

void PairGranHookeHistory::read_restart_settings(FILE *fp) 
{
  if (comm->me == 0) {
    fread(&dampflag,sizeof(int),1,fp);
    fread(&cohesionflag,sizeof(int),1,fp);
  }
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
  MPI_Bcast(&cohesionflag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::reset_dt()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
  returns pointers to properties - for wall/gran 
------------------------------------------------------------------------- */

MechParamGran* PairGranHookeHistory::getMatProp()
{
      if (mpg==NULL) error->all("material properties not defined");
      return mpg;
}
