/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_pour.h"
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

using namespace LAMMPS_NS;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPour::FixPour(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal fix pour command");

  time_depend = 1;

  if (!atom->radius_flag || !atom->rmass_flag)
    error->all("Fix pour requires atom attributes radius, rmass");

  // required args

  ninsert = atoi(arg[3]);
  ntype = atoi(arg[4]);
  seed = atoi(arg[5]);

  if (seed <= 0) error->all("Illegal fix pour command");

  PI = 4.0*atan(1.0);

  // option defaults

  int iregion = -1;
  radius_lo = radius_hi = 0.5;
  density_lo = density_hi = 1.0;
  volfrac = 0.25;
  maxattempt = 50;
  rate = 0.0;
  vxlo = vxhi = vylo = vyhi = vy = vz = 0.0;

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all("Fix pour region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diam") == 0) {
        
      if (iarg+4 > narg) error->all("Illegal fix pour command");
      if(strcmp("uniform",arg[iarg+1])==0) radius_ran_style = RAN_STYLE_UNIFORM;
      else if(strcmp("gaussian",arg[iarg+1])==0) {error->all("Random style gaussian deactivated for radius in fix pour."); radius_ran_style = RAN_STYLE_GAUSSIAN;}
      else if(strcmp("lognormal",arg[iarg+1])==0) radius_ran_style = RAN_STYLE_LOGNORMAL;
      else if(strcmp("discrete",arg[iarg+1])==0) radius_ran_style = RAN_STYLE_DISCRETE;
      else error->all("Illegal fix pour command: Unknown random style");
      radius_lo = 0.5 * atof(arg[iarg+2]);
      radius_hi = 0.5 * atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"dens") == 0) {
        
      if (iarg+4 > narg) error->all("Illegal fix pour command");
      if(strcmp("uniform",arg[iarg+1])==0) density_ran_style = RAN_STYLE_UNIFORM;
      else if(strcmp("gaussian",arg[iarg+1])==0) density_ran_style = RAN_STYLE_GAUSSIAN;
      else if(strcmp("lognormal",arg[iarg+1])==0) density_ran_style = RAN_STYLE_LOGNORMAL;
      else if(strcmp("discrete",arg[iarg+1])==0) density_ran_style = RAN_STYLE_DISCRETE;
      else error->all("Illegal fix pour command: Unknown random style");
      density_lo = atof(arg[iarg+2]);
      density_hi = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"vol") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix pour command");
      volfrac = atof(arg[iarg+1]);
      maxattempt = atoi(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix pour command");
      rate = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (domain->dimension == 3) {
          
	if (iarg+7 > narg) error->all("Illegal fix pour command");
	if(strcmp("uniform",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_UNIFORM;
    else if(strcmp("gaussian",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_GAUSSIAN;
    else if(strcmp("lognormal",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_LOGNORMAL;
    else if(strcmp("discrete",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_DISCRETE;
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
	if(strcmp("uniform",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_UNIFORM;
    else if(strcmp("gaussian",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_GAUSSIAN;
    else if(strcmp("lognormal",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_LOGNORMAL;
    else if(strcmp("discrete",arg[iarg+1])==0) vel_ran_style = RAN_STYLE_DISCRETE;
    else error->all("Illegal fix pour command: Unknown random style");
	vxlo = atof(arg[iarg+2]);
	vxhi = atof(arg[iarg+3]);
	vy = atof(arg[iarg+4]);
	vz = 0.0;
	iarg += 5;
      }
    }
    else if (strcmp(arg[iarg],"template") == 0) iarg+=2; 
    else error->all("Illegal fix pour command");
  }

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

  if(strcmp(this->style,"pour")==0)  
  {
      calc_nfreq();

      force_reneighbor = 1;
      next_reneighbor = update->ntimestep + 1;
      nfirst = next_reneighbor;
      ninserted = 0;

      calc_nper();
  }
}

/* ---------------------------------------------------------------------- */

void FixPour::calc_nfreq()
{
  // nfreq = timesteps between insertions
  // should be time for a particle to fall from top of insertion region
  //   to bottom, taking into account that the region may be moving
  // 1st insertion on next timestep

  double v_relative,delta;
  //double g = 1.; 
  check_gravity(); 
  double g = grav; 

  //fprintf(screen,"grav=%f\n",grav);
  grav=-grav;
 if (domain->dimension == 3) {
   
   v_relative = - sqrt(vz*vz + 2.0*g*(shift_randompos(expectancy(radius_hi,radius_lo,radius_ran_style)))) - rate;
    
   delta = fmax(shift_randompos(expectancy(radius_hi,radius_lo,radius_ran_style)), zhi - zlo -shift_randompos(expectancy(radius_hi,radius_lo,radius_ran_style))) ;

 } else {
   v_relative = vy - rate;
   delta = yhi - ylo;
 }
 double t =   (-v_relative - sqrt(v_relative*v_relative - 2.0*grav*delta)) / grav;
 nfreq = static_cast<int> (t/update->dt + 0.5);

}

/* ---------------------------------------------------------------------- */

void FixPour::calc_nper()
{
  // nper = # to insert each time
  // depends on specified volume fraction
  // volume = volume of insertion region
  // volume_one = average volume of inserted particle
  // in 3d, insure dy >= 1, for quasi-2d simulations

  double volume,volume_one;
  if (domain->dimension == 3) {
    if (region_style == 1) {
      double dy = yhi - ylo;
      if (dy < 1.0) dy = 1.0;
      volume = (xhi-xlo) * dy * (zhi-zlo);
    } else volume = PI*rc*rc * (zhi-zlo);
    volume_one = volume_expectancy(radius_lo,radius_hi,radius_ran_style,3);
  }
  else //2D
  {
    volume = (xhi-xlo) * (yhi-ylo);
    volume_one = volume_expectancy(radius_lo,radius_hi,radius_ran_style,2);
  }

  nper = static_cast<int> (volfrac*volume/volume_one);
  if (nper==0)error->all("Insertion rate too low");

  int nfinal = update->ntimestep + 1 + (ninsert-1)/nper * nfreq;

  // print stats

  if (me == 0) {
    if (screen)
      fprintf(screen,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
    if (logfile)
      fprintf(logfile,
	      "Particle insertion: %d every %d steps, %d by step %d\n",
	      nper,nfreq,ninsert,nfinal);
  }
}

/* ---------------------------------------------------------------------- */

FixPour::~FixPour()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

int FixPour::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPour::check_gravity() 
{
  // insure gravity fix exists
  // for 3d must point in -z, for 2d must point in -y
  // else insertion cannot work

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"gravity") == 0) break;
  if (ifix == modify->nfix)
    error->all("Must use fix gravity before using fix pour");

  double xgrav = ((FixGravity *) modify->fix[ifix])->xgrav*((FixGravity *) modify->fix[ifix])->magnitude;
  double ygrav = ((FixGravity *) modify->fix[ifix])->ygrav*((FixGravity *) modify->fix[ifix])->magnitude;
  double zgrav = ((FixGravity *) modify->fix[ifix])->zgrav*((FixGravity *) modify->fix[ifix])->magnitude;

  if (domain->dimension == 3) {
    if (fabs(xgrav) > EPSILON || fabs(ygrav) > EPSILON ||
	zgrav > EPSILON)
      error->all("Gravity must point in -z to use with fix pour in 3d");
    grav=-zgrav; 
  } else {
    if (fabs(xgrav) > EPSILON || ygrav > EPSILON ||
	fabs(zgrav) > EPSILON)
      error->all("Gravity must point in -y to use with fix pour in 2d");
    grav=-ygrav; 
  }
}

/* ---------------------------------------------------------------------- */

void FixPour::init()
{
  if (domain->triclinic) error->all("Cannot use fix pour with triclinic box");

  check_gravity();

  init_substyle();

}

/* ---------------------------------------------------------------------- */

double FixPour::max_rad(int type)
{
    if (type!=ntype) return 0.;

    if (radius_ran_style==RAN_STYLE_UNIFORM)
    {
        return radius_hi;
    }
    else if (radius_ran_style==RAN_STYLE_GAUSSIAN)
    {
        return radius_lo+3.*radius_hi;
    }
    else if (radius_ran_style==RAN_STYLE_LOGNORMAL)
    {
        error->all("Lognormal expectancy not implemented yet");
    }
    else if (radius_ran_style==RAN_STYLE_DISCRETE)
    {
        error->all("Discrete expectancy not implemented yet");
    }
    return 0.;
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixPour::pre_exchange()
{
  int i;

  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // nnew = # to insert this timestep

  int nnew = nper;
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
    memory->create_2d_double_array(ncount,5,"fix_pour:xmine");  
  double **xnear =
    memory->create_2d_double_array(nprevious+nnew*particles_per_insertion(),5,"fix_pour:xnear");  
  int nnear = nprevious;

  // setup for allgatherv

  int n = 5*ncount; 
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

  MPI_Allgatherv(ptr,5*ncount,MPI_DOUBLE, 
		 xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  // insert new atoms into xnear list, one by one
  // check against all nearby atoms and previously inserted ones
  // if there is an overlap then try again at other z (3d) or y (2d) coord 
  // else insert by adding to xnear list
  // max = maximum # of insertion attempts for all particles
  // h = height, biased to give uniform distribution in time of insertion

  int success;
  double coord[3],radtmp,rn,h;

  int attempt = 0;
  int max = nnew * maxattempt;
  int ntotal = nprevious+nnew;

  while (nnear < ntotal) {
    radtmp = rand_pour(radius_lo,radius_hi,radius_ran_style); 

    success = 0;
    while (attempt < max) {
      
      rn = random->uniform();
      h = (hi_current-shift_randompos(radtmp)) - rn * (hi_current-lo_current-2.*shift_randompos(radtmp));
      attempt++;
      xyz_random(h,coord,radtmp);
      for (i = 0; i < nnear; i++) {
        if(overlaps_xnear_i(coord,radtmp,xnear,i)) break; 
      }
      if (i == nnear) {
	    success = 1;
	    break;
      }
    }
    if (success) {
        nnear=insert_in_xnear(xnear,nnear,coord,radtmp); 
    } else break;
  }

  // warn if not all insertions were performed

  ninserted += nnear-nprevious;
  if (nnear - nprevious < nnew && me == 0)
    error->warning("Less insertions than requested");

  // check if new atom is in my sub-box or above it if I'm highest proc
  // if so, add to my list via create_atom()
  // initialize info about the atom
  // type, diameter, density set from fix parameters
  // group mask set to "all" plus fix group
  // z velocity set to what velocity would be if particle
  //   had fallen from top of insertion region
  // this gives continuous stream of atoms
  // set npartner for new atom to 0 (assume not touching any others)

  AtomVec *avec = atom->avec;
  int j,m,flag;
  double denstmp,vxtmp,vytmp,vztmp;
  double g = grav; 
  //double g = 1.; //originally
  double *sublo = domain->sublo;
  double *subhi = domain->subhi;
  int b_id;  

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  for (i = nprevious; i < nnear; i++) {
    coord[0] = xnear[i][0];
    coord[1] = xnear[i][1];
    coord[2] = xnear[i][2];
    radtmp = xnear[i][3];
    b_id = static_cast<int>(xnear[i][4]);
    denstmp = rand_pour(density_lo,density_hi,density_ran_style)*density_scaling();

    calc_insert_velocities(i,g,xnear,vxtmp,vytmp,vztmp);

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
	  set_body_props(i,m,coord,denstmp,radtmp,vxtmp,vytmp,vztmp,b_id);
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
   insert particle in xnear list
------------------------------------------------------------------------- */

inline int FixPour::insert_in_xnear(double **xnear,int nnear,double *coord,double radtmp)
{
      xnear[nnear][0] = coord[0];
      xnear[nnear][1] = coord[1];
      xnear[nnear][2] = coord[2];
      xnear[nnear][3] = radtmp;
      xnear[nnear][4] = 0.;
      nnear++;
      return nnear;
}

/* ----------------------------------------------------------------------
   check if particle inserted at coord would overlap particle xnear[i]
------------------------------------------------------------------------- */

inline bool FixPour::overlaps_xnear_i(double *coord,double radtmp,double **xnear,int i)
{
    double delx,dely,delz,rsq,radsum;
    delx = coord[0] - xnear[i][0];
	dely = coord[1] - xnear[i][1];
	delz = coord[2] - xnear[i][2];
	rsq = delx*delx + dely*dely + delz*delz;
	radsum = radtmp + xnear[i][3];
	if (rsq <= radsum*radsum) return true;
	else return false;
}

/* ----------------------------------------------------------------------
   how many particles are added to the xnear list in insert_in_xnear()
------------------------------------------------------------------------- */

inline int FixPour::particles_per_insertion()
{
    return 1;
}

/* ----------------------------------------------------------------------
   calculate inlet velocity of particle
------------------------------------------------------------------------- */

inline void FixPour::calc_insert_velocities(int i,double g,double **xnear,double &vxtmp,double &vytmp,double &vztmp)
{
   if (domain->dimension == 3) {
      vxtmp = rand_pour(vxlo,vxhi,vel_ran_style); 
      vytmp = rand_pour(vylo,vyhi,vel_ran_style); 
      //vztmp = vz - sqrt(2.0*g*(hi_current-xnear[i][2]));
      vztmp = - sqrt(vz*vz + 2.0*g*(hi_current-xnear[i][2])); 
    } else {
      vxtmp = rand_pour(vxlo,vxhi,vel_ran_style); 
      vytmp = vy - sqrt(2.0*g*(hi_current-xnear[i][1]));  
      vztmp = 0.0;
    }
}

/* ----------------------------------------------------------------------
   check if particle i could overlap with a particle inserted into region
   return 1 if yes, 0 if no
   use maximum diameter for inserted particle
------------------------------------------------------------------------- */

int FixPour::overlap(int i)
{
  double delta = radius_hi + atom->radius[i];
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

double FixPour::rand_pour(double param1, double param2, int style)
{
    if (style==RAN_STYLE_UNIFORM)
    {
        return param1 + random->uniform() * (param2-param1);
    }
    else if (style==RAN_STYLE_GAUSSIAN)
    {
        return param2 * random->gaussian() + param1;
    }
    else if (style==RAN_STYLE_LOGNORMAL)
    {
        error->all("Lognormal random not implemented yet");
    }
    else if (style==RAN_STYLE_DISCRETE)
    {
        error->all("Discrete random not implemented yet");
    }
}

/* ---------------------------------------------------------------------- */

//uniform:  param 1 = low  param 2 = high
//gaussian: param 1 = mu   param 2 = sigma
double FixPour::expectancy(double param1, double param2, int style)
{
    if (style==RAN_STYLE_UNIFORM)
    {
        return (param2+param1)/2.;
    }
    else if (style==RAN_STYLE_GAUSSIAN)
    {
        return param1;
    }
    else if (style==RAN_STYLE_LOGNORMAL)
    {
        error->all("Lognormal expectancy not implemented yet");
    }
    else if (style==RAN_STYLE_DISCRETE)
    {
        error->all("Discrete expectancy not implemented yet");
    }
    return 0.;
}

/* ---------------------------------------------------------------------- */

//uniform:  param 1 = low  param 2 = high
//gaussian: param 1 = mu   param 2 = sigma
double FixPour::volume_expectancy(double param_1, double param_2, int style,int dim)
{
    if (style==RAN_STYLE_UNIFORM)
    {
        if(dim==3) //3D
        {
            if(param_2/param_1<0.9999) error->all("fix pour: radius_hi must be larger than radius_lo");
            if(param_2/param_1<1.0002) return (4.*PI/3.* param_2*param_2*param_2);
            else return (PI/3.0 * (pow(param_2,4.)-pow(param_1,4.))/(param_2-param_1));
        }
        else //2D
        {
            if(param_2/param_1<0.9999) error->all("fix pour: radius_hi must be larger than radius_lo");
            if(param_2/param_1<1.0002) return (4.*PI* param_2*param_2);
            return (PI/3.0 * (pow(param_2,3.)-pow(param_1,3.))/(param_2-param_1));
        }
    }
    else if (style==RAN_STYLE_GAUSSIAN)
    {
        if(dim==3) //3D
            return (4.*PI/3.*param_1*(param_1*param_1+3.*param_2*param_2));
        else //2D
            return ( PI * (param_1*param_1 + param_2*param_2));
    }
    else if (style==RAN_STYLE_LOGNORMAL)
        error->all("Lognormal expectancy not implemented yet");
    else if (style==RAN_STYLE_DISCRETE)
        error->all("Discrete expectancy not implemented yet");
    return 0.;
}

/* ---------------------------------------------------------------------- */

inline double FixPour::density_scaling()
{
    return 1.;
}

/* ---------------------------------------------------------------------- */

inline double FixPour::shift_randompos(double rt)
{
    return rt;
}

/* ---------------------------------------------------------------------- */

void FixPour::xyz_random(double h, double *coord,double rt) 
{
  if (domain->dimension == 3) {
    if (region_style == 1) {
      coord[0] = xlo+shift_randompos(rt) + random->uniform() * (xhi-xlo-2.*shift_randompos(rt));  
      coord[1] = ylo+shift_randompos(rt) + random->uniform() * (yhi-ylo-2.*shift_randompos(rt));  
      coord[2] = h;
    } else {
      double r1,r2;
      while (1) {
	r1 = random->uniform() - 0.5;
	r2 = random->uniform() - 0.5;
	if (r1*r1 + r2*r2 < 0.25) break;
      }
      coord[0] = xc + 2.0*r1*(rc-shift_randompos(rt)); 
      coord[1] = yc + 2.0*r2*(rc-shift_randompos(rt)); 
      coord[2] = h;
    }
  } else {
    coord[0] = xlo+shift_randompos(rt) + random->uniform() * (xhi-xlo-2.*shift_randompos(rt)); 
    coord[1] = h;
    coord[2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPour::reset_dt()
{
  error->all("Cannot change timestep with fix pour");
}
