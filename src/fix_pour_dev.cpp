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
#include "fix_pour_dev.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "domain.h"
#include "region.h"
#include "region_block.h"
#include "region_cylinder.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "particleToInsert.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixPourDev::FixPourDev(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all("Illegal fix pour command, not enough arguments");

  time_depend = 1;
  restart_global = 1; 

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all("Fix pour requires atom attributes radius, rmass");

  ninsert = 0;
  masstotal = 0.;
  mass_ins = 0.;

  nBody=0; 

  check_ol_flag = 1;

  // required args

  //one of the two criteria is reduced 
  iarg=3;

  if(strcmp(this->style,"pour/dev/packing")!=0)
  {
      if(strcmp(arg[iarg],"nparticles")==0) ninsert = atoi(arg[iarg+1]);
      else if(strcmp(arg[iarg],"mass")==0) masstotal = atof(arg[iarg+1]);
      else error->all("Fix pour has to specify either target number or target mass of particles to insert");
      iarg+=2;
  }

  seed = atoi(arg[iarg++]);
  if (seed <= 0) error->all("Illegal fix pour command");
  PI = 4.0*atan(1.0);

  //get template
  if(strcmp(arg[iarg++],"distributiontemplate")!=0)  error->all("Fix pour command requires you to define a distributiontemplate");
  int ifix=modify->find_fix(arg[iarg++]);
  if(ifix<0) error->all("Fix pour command requires you to define a valid ID for a fix of type particledistribution/discrete");
  if(strcmp(modify->fix[ifix]->style,"particledistribution/discrete")!=0) error->all("Fix pour command requires you to define a valid ID for a fix of type particledistribution/discrete");
  fpdd=static_cast<FixParticledistributionDiscrete*>(modify->fix[ifix]);

  //check appropriateness of fix pour
  if(strcmp(this->style,"pour/dev")==0&&fpdd->max_nspheres()>1) error->all("Your discrete particle distribution defines multisphere particles. Thus, you have to use fix pour/multisphere");

  //min/max type to be inserted, need that to check if material properties defined for all materials
  ntype_max=fpdd->max_type();
  ntype_min=fpdd->min_type();

  //calculate ninsert from total mass and template expectancy 
  //ninsert is now enforcing total mass if specified
  if(masstotal==0) masstotal=ninsert*fpdd->mass_expect();
  if(ninsert==0) ninsert = MAX(ninsert,static_cast<int>(masstotal/fpdd->mass_expect()));

  if(ninsert<1. && strcmp(this->style,"pour/dev/packing")!=0) error->all("Can not insert less than one particle");

  // option defaults

  int iregion = -1;
  volfrac = 0.25;
  massflowrate = 0.0; 
  maxattempt = 50;
  rate = 0.0;
  vel_rand_style = RAN_STYLE_UNIFORM_FP;
  vxlo = vxhi = vylo = vyhi = vy = vz = 0.0;
  nRegEx = 0;

  // optional args

  //iarg = 8; 
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all("Fix pour region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix pour command");
      volfrac = atof(arg[iarg+1]);
      maxattempt = atoi(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"massflowrate") == 0) {   
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      massflowrate = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      rate = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
          if (domain->dimension == 3) {
              
            if (iarg+7 > narg) error->all("Illegal fix pour command");
            if(strcmp("uniform",arg[iarg+1])==0) vel_rand_style = RAN_STYLE_UNIFORM_FP;
            else if(strcmp("gaussian",arg[iarg+1])==0) vel_rand_style = RAN_STYLE_GAUSSIAN_FP;
            else error->all("Illegal fix pour command: Unknown random style");
            vxlo = atof(arg[iarg+2]);
            vxhi = atof(arg[iarg+3]);
            vylo = atof(arg[iarg+4]);
            vyhi = atof(arg[iarg+5]);
            vz = atof(arg[iarg+6]);
            if(vz>0.) error->all("fix pour: z-velocity must not not be <0.");
            iarg += 7;
          } else {
          
            if (iarg+5 > narg) error->all("Illegal fix pour command");
            if(strcmp("uniform",arg[iarg+1])==0) vel_rand_style = RAN_STYLE_UNIFORM_FP;
            else if(strcmp("gaussian",arg[iarg+1])==0) vel_rand_style = RAN_STYLE_GAUSSIAN_FP;
            else error->all("Illegal fix pour command: Unknown random style");
            vxlo = atof(arg[iarg+2]);
            vxhi = atof(arg[iarg+3]);
            vy = atof(arg[iarg+4]);
            vz = 0.0;
            iarg += 5;
           }
    } else if (strcmp(arg[iarg],"regionexempts") == 0) { 
      if (iarg+3 > narg) error->all("Illegal fix pour command");
      nRegEx=atoi(arg[iarg+1]);
      if(nRegEx<1) error->all("Fix pour: Number of region exempts must be > 0");
      if (iarg+2+nRegEx > narg) error->all("Illegal fix pour command");
      regExList=new Region*[nRegEx];
      for (int i=0;i<nRegEx;i++)
      {
          int eregion = domain->find_region(arg[iarg+2+i]);
          if (eregion == -1) error->all("Fix pour region ID does not exist");
          regExList[i]=domain->regions[eregion];
      }
      iarg += 2+nRegEx;
    } else if (strcmp(arg[iarg],"overlapcheck") == 0) { 
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      if(strcmp(arg[iarg+1],"yes")==0) check_ol_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) check_ol_flag = 0;
      else error->all("Illegal fix pour command");
      iarg += 2;;
    }
    else error->all("Illegal fix pour command");
  }

  if(massflowrate!=0.&&comm->me==0&&screen)  fprintf(screen, "INFO: You are specifying a massflow rate. This may result in a lower volume fraction than specified.\n");
  if(massflowrate!=0.&&comm->me==0&&logfile) fprintf(logfile,"INFO: You are specifying a massflow rate. This may result in a lower volume fraction than specified.\n");

  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all("Must specify a region in fix pour");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all("Fix pour region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all("Fix pour region cannot be dynamic");

  if (strcmp(domain->regions[iregion]->style,"block") == 0) {
    region_style = 1;
    xlo = ((RegBlock *) domain->regions[iregion])->xlo;
    xhi = ((RegBlock *) domain->regions[iregion])->xhi;
    ylo = ((RegBlock *) domain->regions[iregion])->ylo;
    yhi = ((RegBlock *) domain->regions[iregion])->yhi;
    zlo = ((RegBlock *) domain->regions[iregion])->zlo;
    zhi = ((RegBlock *) domain->regions[iregion])->zhi;
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
	ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
	zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all("Insertion region extends outside simulation box");
  } else if (strcmp(domain->regions[iregion]->style,"cylinder") == 0) {
    region_style = 2;
    char axis = ((RegCylinder *) domain->regions[iregion])->axis;
    xc = ((RegCylinder *) domain->regions[iregion])->c1;
    yc = ((RegCylinder *) domain->regions[iregion])->c2;
    rc = ((RegCylinder *) domain->regions[iregion])->radius;
    zlo = ((RegCylinder *) domain->regions[iregion])->lo;
    zhi = ((RegCylinder *) domain->regions[iregion])->hi;
    if (axis != 'z')
      error->all("Must use a z-axis cylinder with fix pour");
    if (xc-rc < domain->boxlo[0] || xc+rc > domain->boxhi[0] ||
	yc-rc < domain->boxlo[1] || yc+rc > domain->boxhi[1] ||
	zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all("Insertion region extends outside simulation box");
  } else error->all("Must use a block or cylinder region with fix pour");

  if (region_style == 2 && domain->dimension == 2)
    error->all("Must use a block region with fix pour for 2d simulations");

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);

  // allgather arrays

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // 1st insertion on next timestep

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;

  // nper = # to insert each time
  // depends on specified volume fraction
  // volume = volume of insertion region
  // volume_one = volume of inserted particle (with max possible radius)
  // in 3d, insure dy >= 1, for quasi-2d simulations

  double volume,volume_one;
  if (domain->dimension == 3) {
    if (region_style == 1) {
      double dy = yhi - ylo;
      
      volume = (xhi-xlo) * dy * (zhi-zlo);
    } else volume = PI*rc*rc * (zhi-zlo);
    volume_one = fpdd->vol_expect(); 
  } else {
    volume = (xhi-xlo) * (yhi-ylo);
    volume_one = fpdd->vol_expect(); 
  }

  nper = volfrac*volume/volume_one;       
  nfinal = static_cast<int>(update->ntimestep + 1 + (ninsert-1)/nper * nfreq); 

  //return here for doing packing
  if(strcmp(this->style,"pour/dev/packing")==0) return;

  // grav = gravity in distance/time^2 units
  // assume grav = -magnitude at this point, enforce in init()

  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
  if (ifix == modify->nfix)
    error->all("No fix gravity defined for fix pour");
  grav = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;

  // nfreq = timesteps between insertions
  // should be time for a particle to fall from top of insertion region
  //   to bottom, taking into account that the region may be moving
  // set these 2 eqs equal to each other, solve for smallest positive t
  //   x = zhi + vz*t + 1/2 grav t^2
  //   x = zlo + rate*t
  //   gives t = [-(vz-rate) - sqrt((vz-rate)^2 - 2*grav*(zhi-zlo))] / grav
  //   where zhi-zlo > 0, grav < 0, and vz & rate can be either > 0 or < 0

  //account for the fact that particles have starting velocities v0_add, due to the
  //fact that insertion starts a bit lower
  //use diameter of smallest particle class as reference

  double v_relative,delta;
  double v0_add=-vz-sqrt(vz*vz-2.0*grav*(fpdd->templates[fpdd->distorder[fpdd->ntemplates-1]]->pti->r_bound));
  
  if (domain->dimension == 3) {
    v_relative = vz - rate;
    delta = zhi - zlo;
  } else {
    v_relative = vy - rate;
    delta = yhi - ylo;
  }
  v_relative+=v0_add;
  double t =
    (-v_relative - sqrt(v_relative*v_relative - 2.0*grav*delta)) / grav;
  nfreq = static_cast<int> (t/update->dt + 0.5);

  // nper2 = # to insert each time 
  // only choose if massflowrate is criterion; depends on specified mass flow rate and mass expectancy
  if(massflowrate>0.)
  {
      double nper2 = nfreq*update->dt*massflowrate/fpdd->mass_expect();
      if(nper2<1.) error->all("Would insert less than 1 particle per insertion step, please increase mass-flow rate or make insertion volume larger");
      int nfinal2 = static_cast<int>(update->ntimestep + 1 + (ninsert-1)/nper2 * nfreq);

      //choose min of nper1 and nper2
      nper = MIN (nper,nper2);
      if(nper==nper2) nfinal=nfinal2;
      //if volume fraction is the restrictive criterion, warn because defined mass flow rate will not be reached
      else if(me==0) error->all("Fix pour/dev can not achieve desired mass flow rate, choose higher volume fraction or higher insertion velocity");
  }

  // print stats

  if (me == 0) {
    if (screen)
      if(particles_per_insertion()==1) fprintf(screen,
	      "Particle insertion: %f every %d steps, %d within %d steps\n",
	      nper,nfreq,ninsert,nfinal-update->ntimestep);
	  else fprintf(screen,
	      "Rigid body insertion: %f every %d steps, %d within %d steps\nEach rigid body consists of %d particles\n",
	      nper,nfreq,ninsert,nfinal-update->ntimestep,particles_per_insertion());
    if (logfile)
      if(particles_per_insertion()==1) fprintf(logfile,
	      "Particle insertion: %f every %d steps, %d within %d steps\n",
	      nper,nfreq,ninsert,nfinal-update->ntimestep);
      else fprintf(logfile,
	      "Rigid body insertion: %f every %d steps, %d within %d steps\nEach rigid body consists of %d particles\n",
	      nper,nfreq,ninsert,nfinal-update->ntimestep,particles_per_insertion());
  }

}

/* ---------------------------------------------------------------------- */

FixPourDev::~FixPourDev()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;

  if(nRegEx>0) delete []regExList;
}

/* ---------------------------------------------------------------------- */

int FixPourDev::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPourDev::init()
{
  if (domain->triclinic) error->all("Cannot use fix pour with triclinic box");

  // insure gravity fix exists
  // for 3d must point in -z, for 2d must point in -y
  // else insertion cannot work

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
  if (ifix == modify->nfix)
    error->all("No fix gravity defined for fix pour");

  double xgrav = ((FixGravity *) modify->fix[ifix])->xgrav;
  double ygrav = ((FixGravity *) modify->fix[ifix])->ygrav;
  double zgrav = ((FixGravity *) modify->fix[ifix])->zgrav;

  if (domain->dimension == 3) {
    if (fabs(xgrav) > EPSILON || fabs(ygrav) > EPSILON ||
	fabs(zgrav+1.0) > EPSILON)
      error->all("Gravity must point in -z to use with fix pour in 3d");
  } else {
    if (fabs(xgrav) > EPSILON || fabs(ygrav+1.0) > EPSILON ||
	fabs(zgrav) > EPSILON)
      error->all("Gravity must point in -y to use with fix pour in 2d");
  }

  double gnew = - ((FixGravity *) modify->fix[ifix])->magnitude * force->ftm2v;
  if (gnew != grav)
    error->all("Gravity changed since fix pour was created");
  init_substyle();
}

/* ---------------------------------------------------------------------- */

double FixPourDev::max_rad(int type)
{
    return fpdd->max_rad(type);
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixPourDev::pre_exchange()
{

  int i;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // nnew = # to insert this timestep

  int nnew = static_cast<int>(nper+random->uniform());

  //init number of bodies to be inserted
  nnew = fpdd->random_init(nnew);

  if (ninserted + nnew > ninsert) nnew = ninsert - ninserted;

  // lo/hi current = z (or y) bounds of insertion region this timestep

  if (domain->dimension == 3) {
    lo_current = zlo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = zhi + (update->ntimestep - nfirst) * update->dt * rate;
  } else {
    lo_current = ylo + (update->ntimestep - nfirst) * update->dt * rate;
    hi_current = yhi + (update->ntimestep - nfirst) * update->dt * rate;
  }

  // ncount = # of my atoms that overlap the insertion region
  // nprevious = total of ncount across all procs

  int ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) ncount++;

  int nprevious;
  MPI_Allreduce(&ncount,&nprevious,1,MPI_INT,MPI_SUM,world);

  // xmine is for my atoms
  // xnear is for atoms from all procs + atoms to be inserted

  double **xmine =
    memory->create_2d_double_array(ncount,7,"fix_pour_dev:xmine");
  double **xnear =
    memory->create_2d_double_array(nprevious+nnew*particles_per_insertion(),7,"fix_pour_dev:xnear");
  int nnear = nprevious;

  // setup for allgatherv

  int n = 7*ncount;
  MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // load up xmine array

  double **x = atom->x;
  double *radius = atom->radius;

  ncount = 0;
  for (i = 0; i < atom->nlocal; i++)
    if (overlap(i)) {
      xmine[ncount][0] = x[i][0];
      xmine[ncount][1] = x[i][1];
      xmine[ncount][2] = x[i][2];
      xmine[ncount][3] = radius[i];
      ncount++;
    }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = NULL;
  if (ncount) ptr = xmine[0];
  MPI_Allgatherv(ptr,7*ncount,MPI_DOUBLE,
		 xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  // insert new atoms into xnear list, one by one
  // check against all nearby atoms and previously inserted ones
  // if there is an overlap then try again at same z (3d) or y (2d) coord
  // else insert by adding to xnear list
  // max = maximum # of insertion attempts for all particles
  // h = height, biased to give uniform distribution in time of insertion

  int success;
  int isInExempt;
  double coord[3],radtmp,delx,dely,delz,rsq,radsum,rn,h,tmp,vzrel;

  int attempt = 0;
  int max = nnew * maxattempt;

  int ntotal = nprevious+nnew;

  vzrel=vz-rate;

  if(hi_current-lo_current<3.*shift_randompos()) error->all("Fix pour insertion volume is not tall enough for largest insertion, make it at least 3 particle radii tall.");

  while (nnear < ntotal) {

    //randomize particle properties
    fpdd->randomize();
    isInExempt = 0; 

    do
    {
        random_insert_height(h,tmp,vzrel);
        
    }
    while(h<(lo_current+shift_randompos())||h>(hi_current-shift_randompos()));

    success = 0;
    while (attempt < max) {
          attempt++;
          xyz_random(h,coord);
          if(isInExempt = isInExemptRegion(coord) == 1) break; 

          if(!check_ol_flag)
          {
              success = 1;
              break;
          }

          for (i = 0; i < nnear; i++) {
            if(overlaps_xnear_i(coord,xnear,i)) break; 
          }
          if (i == nnear) {
            success = 1;
            break;
          }
    }

    if (success) {
      nnear=insert_in_xnear(xnear,nnear,coord); 
      //if more that 1 particle per body added, add number of particles -1
      ntotal+=(fpdd->pti->nspheres-1);
      ninserted ++;
    } else if(attempt == max)break; 
  }

  // warn if not all insertions were performed

  if (nnear - nprevious < nnew && me == 0)
    error->warning("Less insertions than requested, the particle size distribution you wish may not be pictured");

  //show comparison of ideal massflowrate and calculated one 
  if(massflowrate>0. && me==0)
  {
      double ts=static_cast<double>(MIN(update->ntimestep-nfirst+nfreq,nfinal-nfirst+nfreq));
      if(screen)  fprintf(screen, "Mass inserted: %f (expectancy value: %f)\n",mass_ins,fmin(masstotal,update->dt*ts*massflowrate));
      if(logfile) fprintf(logfile,"Mass inserted: %f (expectancy value: %f)\n",mass_ins,fmin(masstotal,update->dt*ts*massflowrate));
  }

  // check if new atom is in my sub-box or above it if I'm highest proc
  // if so, add to my list via create_atom()
  // initialize info about the atom
  // type, diameter, density set from fix parameters
  // group mask set to "all" plus fix group
  // z velocity set to what velocity would be if particle
  //   had fallen from top of insertion region
  //   this gives continuous stream of atoms
  //   solution for v from these 2 eqs, after eliminate t:
  //     v = vz + grav*t
  //     coord[2] = hi_current + vz*t + 1/2 grav t^2
  // set npartner for new atom to 0 (assume not touching any others)

  AtomVec *avec = atom->avec;
  int j,m,flag;
  double denstmp,vxtmp,vytmp,vztmp;
  int ntype,mol_id;
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  for (i = nprevious; i < nnear; i++) {
    coord[0] = xnear[i][0];
    coord[1] = xnear[i][1];
    coord[2] = xnear[i][2];
    radtmp = xnear[i][3];
    denstmp = xnear[i][4];
    mol_id = static_cast<int>(xnear[i][5]);
    ntype = static_cast<int>(xnear[i][6]);
    calc_insert_velocities(i,xnear,vxtmp,vytmp,vztmp);

    flag = 0;
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	coord[1] >= sublo[1] && coord[1] < subhi[1] &&
	coord[2] >= sublo[2] && coord[2] < subhi[2]) flag = 1;
    else if (domain->dimension == 3 && coord[2] >= domain->boxhi[2] &&
	     comm->myloc[2] == comm->procgrid[2]-1 &&
	     coord[0] >= sublo[0] && coord[0] < subhi[0] &&
	     coord[1] >= sublo[1] && coord[1] < subhi[1]) flag = 1;
    else if (domain->dimension == 2 && coord[1] >= domain->boxhi[1] &&
	     comm->myloc[1] == comm->procgrid[1]-1 &&
	     coord[0] >= sublo[0] && coord[0] < subhi[0]) flag = 1;

    if (flag) {
      
      avec->create_atom(ntype,coord);
      m = atom->nlocal - 1;
      atom->type[m] = ntype;
      atom->radius[m] = radtmp;
      atom->density[m] = denstmp;
      atom->rmass[m] = 4.0*PI/3.0 * radtmp*radtmp*radtmp * denstmp;
      atom->mask[m] = 1 | groupbit;
      atom->v[m][0] = vxtmp;
      atom->v[m][1] = vytmp;
      atom->v[m][2] = vztmp;

      for (j = 0; j < nfix; j++)
	    if (fix[j]->create_attribute) fix[j]->set_arrays(m);
    }
	if(give_mol_id()&&mol_id>=0) { 
          set_body_props(flag,m,vxtmp,vytmp,vztmp,mol_id); 
    }
  }

  // set tag # of new particles beyond all previous atoms
  // reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  if (atom->tag_enable) {
    atom->tag_extend();
    atom->natoms += nnear - nprevious;
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  finalize_insertion();

  // free local memory

  memory->destroy_2d_double_array(xmine);
  memory->destroy_2d_double_array(xnear);

  // next timestep to insert

  if (ninserted < ninsert) next_reneighbor += nfreq;
  else next_reneighbor = 0;
}

/* ----------------------------------------------------------------------
   check if particle i could overlap with a particle inserted into region
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int FixPourDev::overlap(int i)
{
  double delta = atom->radius[i]; 
  double **x = atom->x;

  if (domain->dimension == 3) {
    if (region_style == 1) {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < ylo-delta || x[i][1] > yhi+delta ||
	  x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
    } else {
      if (x[i][2] < lo_current-delta || x[i][2] > hi_current+delta) return 0;
      double delx = x[i][0] - xc;
      double dely = x[i][1] - yc;
      double rsq = delx*delx + dely*dely;
      double r = rc + delta;
      if (rsq > r*r) return 0;
    }
  } else {
      if (x[i][0] < xlo-delta || x[i][0] > xhi+delta ||
	  x[i][1] < lo_current-delta || x[i][1] > hi_current+delta) return 0;
  }

  return 1;
}

/* ---------------------------------------------------------------------- */

void FixPourDev::xyz_random(double h, double *coord)
{
  if (domain->dimension == 3) {
    if (region_style == 1) {
      coord[0] = xlo+shift_randompos() + random->uniform() * (xhi-xlo-2.*shift_randompos());  
      coord[1] = ylo+shift_randompos() + random->uniform() * (yhi-ylo-2.*shift_randompos());  
      coord[2] = h;
    } else {
      double r1,r2;
      while (1) {
	r1 = random->uniform() - 0.5;
	r2 = random->uniform() - 0.5;
	if (r1*r1 + r2*r2 < 0.25) break;
      }
      coord[0] = xc + 2.0*r1*(rc-shift_randompos()); 
      coord[1] = yc + 2.0*r2*(rc-shift_randompos()); 
      coord[2] = h;
    }
  } else {
    coord[0] = xlo+shift_randompos() + random->uniform() * (xhi-xlo-2.*shift_randompos()); 
    coord[1] = h;
    coord[2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

double FixPourDev::shift_randompos()
{
    return fpdd->pti->r_bound;
}

/* ---------------------------------------------------------------------- */

void FixPourDev::reset_dt()
{
  error->all("Cannot change timestep with fix pour");
}

/* ----------------------------------------------------------------------
   test inserted particles versus particles to insert
------------------------------------------------------------------------- */
inline bool FixPourDev::overlaps_xnear_i(double *coord,double **xnear,int i)
{
    double delx,dely,delz,rsq,radsum;
    for(int j=0;j<fpdd->pti->nspheres;j++)
    {
       delx = coord[0] + fpdd->pti->x_ins[j][0] - xnear[i][0];
	   dely = coord[1] + fpdd->pti->x_ins[j][1] - xnear[i][1];
	   delz = coord[2] + fpdd->pti->x_ins[j][2] - xnear[i][2];
	   rsq = delx*delx + dely*dely + delz*delz;
	   radsum = fpdd->pti->radius_ins[j] + xnear[i][3];
	   if (rsq <= radsum*radsum) return true;
    }
	return false;
}

/* ----------------------------------------------------------------------
   insert particle in xnear list
------------------------------------------------------------------------- */

inline int FixPourDev::insert_in_xnear(double **xnear,int nnear,double *coord)
{
    
    for(int j=0;j<fpdd->pti->nspheres;j++)
    {
      xnear[nnear][0] = coord[0] + fpdd->pti->x_ins[j][0];
      xnear[nnear][1] = coord[1] + fpdd->pti->x_ins[j][1];
      xnear[nnear][2] = coord[2] + fpdd->pti->x_ins[j][2];
      xnear[nnear][3] = fpdd->pti->radius_ins[j];
      xnear[nnear][4] = fpdd->pti->density_ins;
      xnear[nnear][5] = static_cast<double>(nBody);
      xnear[nnear][6] = static_cast<double>(fpdd->pti->atom_type);
      
      nnear++;
    }
    nBody++;
    mass_ins+=fpdd->pti->mass_ins;
    return nnear;
}

/* ----------------------------------------------------------------------
   calculate inlet velocity of particle
------------------------------------------------------------------------- */

inline void FixPourDev::calc_insert_velocities(int i,double **xnear,double &vxtmp,double &vytmp,double &vztmp)
{
   if (domain->dimension == 3) {
      vxtmp = rand_pour(vxlo,vxhi,vel_rand_style); 
      vytmp = rand_pour(vylo,vyhi,vel_rand_style); 
      vztmp = - sqrt(vz*vz - 2.0*grav*(hi_current-xnear[i][2])); 
    } else {
      vxtmp = rand_pour(vxlo,vxhi,vel_rand_style); 
      vytmp = vy - sqrt(2.0*grav*(hi_current-xnear[i][1]));  
      vztmp = 0.0;
    }
}

/* ----------------------------------------------------------------------
   random insertion heights of particle
------------------------------------------------------------------------- */

inline void FixPourDev::random_insert_height(double &h,double &tmp,double vzrel)
{
        double rn = random->uniform();
        tmp=rn*(vzrel+sqrt(-2.*grav*(hi_current-lo_current)+vzrel*vzrel))-vzrel;
        h = (tmp*tmp-vzrel*vzrel)/(-2.*grav);
        h = hi_current - h;
}

/* ---------------------------------------------------------------------- */

double FixPourDev::rand_pour(double param1, double param2, int style)
{
    if (style==RAN_STYLE_UNIFORM_FP)
    {
        return param1 + random->uniform() * (param2-param1);
    }
    else if (style==RAN_STYLE_GAUSSIAN_FP)
    {
        return param2 * random->gaussian() + param1;
    }
    else if (style==RAN_STYLE_CONSTANT_FP)
    {
        return param1;
    }
}

/* ---------------------------------------------------------------------- */

inline int FixPourDev::isInExemptRegion(double * coo)
{
    if(nRegEx==0) return 0;

    int inside=0;

    for(int i=0;i<nRegEx;i++)
    {
        if(regExList[i]->inside(coo[0],coo[1],coo[2])) inside++;
    }
    return inside;
}

/* ---------------------------------------------------------------------- */

inline bool FixPourDev::give_mol_id()
{
    if(fpdd->max_nspheres()>1) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

inline int FixPourDev::particles_per_insertion()
{
   //worst case: template with the most spheres is chosen every time
    return fpdd->max_nspheres();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixPourDev::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = static_cast<double>(random->state());
  list[n++] = static_cast<double>(ninserted);
  list[n++] = static_cast<double>(nfirst);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = static_cast<double>(nBody);
  list[n++] = mass_ins;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixPourDev::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  ninserted = static_cast<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);
  nBody = static_cast<int> (list[n++]);
  mass_ins= list[n++];

  random->reset(seed);
}
