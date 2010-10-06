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
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "modify.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "fix_move_tri.h"
#include "fix_meshGran.h"

using namespace LAMMPS_NS;

enum{LINEAR,WIGGLE,RIGGLE,ROTATE,DEFORM,VARIABLE};
enum{EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixMoveTri::FixMoveTri(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if (screen && comm->me==0) fprintf(screen,"Note: A fix move/mesh/gran error or warning may appear as fix move error or warning\n");

  if (narg < 8) error->all("Illegal fix move/mesh/gran command, not enough arguments");
  if (!atom->radius_flag || !atom->rmass_flag)
    error->all("Fix move/mesh/gran requires atom attributes radius, rmass (atom style granular)");

  //the last arg must be a valid skin safety factor 
  skinSafetyFactor=atoi(arg[narg-1]);

  //the arg before must be a valid fix of style fix mesh/gran 
  //if not found, assume that the last arg is the reference to the fix mesh gran
  char* f_id=arg[narg-2];
  int f_i=modify->find_fix(f_id);
  if (f_i==-1)
  {
      f_id=arg[narg-1];
      f_i=modify->find_fix(f_id);
  }

  if (f_i==-1) error->all("Can not find the fix mesh/gran to apply the fix move/mesh/gran to");
  if (strncmp(modify->fix[f_i]->style,"mesh/gran",8)!=0) error->all("Can apply fix move/mesh/gran only to a fix of type fix mesh/gran");

  time_integrate = 0; 
  time_depend = 1;

  // parse args

  int iarg;
  xvarstr = yvarstr = zvarstr = NULL;
  vxvarstr = vyvarstr = vzvarstr = NULL;

  if (strcmp(arg[3],"linear") == 0) {
    if (narg < 7) error->all("Illegal fix move command");
    iarg = 7;
    mstyle = LINEAR;
    if (strcmp(arg[4],"NULL") == 0) vxflag = 0;
    else {
      vxflag = 1;
      vx = atof(arg[4]);
    }
    if (strcmp(arg[5],"NULL") == 0) vyflag = 0;
    else {
      vyflag = 1;
      vy = atof(arg[5]);
    }
    if (strcmp(arg[6],"NULL") == 0) vzflag = 0;
    else {
      vzflag = 1;
      vz = atof(arg[6]);
    }

  } else if (strcmp(arg[3],"wiggle") == 0) {
    if (narg < 8) error->all("Illegal fix move command");
    iarg = 8;
    mstyle = WIGGLE;
    if (strcmp(arg[4],"NULL") == 0) axflag = 0;
    else {
      axflag = 1;
      ax = atof(arg[4]);
    }
    if (strcmp(arg[5],"NULL") == 0) ayflag = 0;
    else {
      ayflag = 1;
      ay = atof(arg[5]);
    }
    if (strcmp(arg[6],"NULL") == 0) azflag = 0;
    else {
      azflag = 1;
      az = atof(arg[6]);
    }
    period = atof(arg[7]);

  } else if (strcmp(arg[3],"riggle") == 0) {  
    if (narg < 12) error->all("Illegal fix move command");
    iarg = 12;
    mstyle = RIGGLE;
    point[0] = atof(arg[4]);
    point[1] = atof(arg[5]);
    point[2] = atof(arg[6]);
    axis[0] = atof(arg[7]);
    axis[1] = atof(arg[8]);
    axis[2] = atof(arg[9]);
    period = atof(arg[10]);
	ampl = atof(arg[11]);

  } else if (strcmp(arg[3],"rotate") == 0) {
    if (narg < 11) error->all("Illegal fix move command");
    iarg = 11;
    mstyle = ROTATE;
    point[0] = atof(arg[4]);
    point[1] = atof(arg[5]);
    point[2] = atof(arg[6]);
    axis[0] = atof(arg[7]);
    axis[1] = atof(arg[8]);
    axis[2] = atof(arg[9]);
    period = atof(arg[10]);

  } else if (strcmp(arg[3],"deform/periodic") == 0) { 
    if (narg < 9) error->all("Illegal fix move command");
    iarg = 9;
    mstyle = DEFORM;
    point[0] = atof(arg[4]);
    point[1] = atof(arg[5]);
    point[2] = atof(arg[6]);
    period = atof(arg[7]);
    ampl = atof(arg[8]);
  } else if (strcmp(arg[3],"variable") == 0) {
    if (narg < 10) error->all("Illegal fix move command");
    iarg = 10;
    mstyle = VARIABLE;
    if (strcmp(arg[4],"NULL") == 0) xvarstr = NULL;
    else if (strstr(arg[4],"v_") == arg[4]) {
      int n = strlen(&arg[4][2]) + 1;
      xvarstr = new char[n];
      strcpy(xvarstr,&arg[4][2]);
    } else error->all("Illegal fix move command");
    if (strcmp(arg[5],"NULL") == 0) yvarstr = NULL;
    else if (strstr(arg[5],"v_") == arg[5]) {
      int n = strlen(&arg[5][2]) + 1;
      yvarstr = new char[n];
      strcpy(yvarstr,&arg[5][2]);
    } else error->all("Illegal fix move command");
    if (strcmp(arg[6],"NULL") == 0) zvarstr = NULL;
    else if (strstr(arg[6],"v_") == arg[6]) {
      int n = strlen(&arg[6][2]) + 1;
      zvarstr = new char[n];
      strcpy(zvarstr,&arg[6][2]);
    } else error->all("Illegal fix move command");
    if (strcmp(arg[7],"NULL") == 0) vxvarstr = NULL;
    else if (strstr(arg[7],"v_") == arg[7]) {
      int n = strlen(&arg[7][2]) + 1;
      vxvarstr = new char[n];
      strcpy(vxvarstr,&arg[7][2]);
    } else error->all("Illegal fix move command");
    if (strcmp(arg[8],"NULL") == 0) vyvarstr = NULL;
    else if (strstr(arg[8],"v_") == arg[8]) {
      int n = strlen(&arg[8][2]) + 1;
      vyvarstr = new char[n];
      strcpy(vyvarstr,&arg[8][2]);
    } else error->all("Illegal fix move command");
    if (strcmp(arg[9],"NULL") == 0) vzvarstr = NULL;
    else if (strstr(arg[9],"v_") == arg[9]) {
      int n = strlen(&arg[9][2]) + 1;
      vzvarstr = new char[n];
      strcpy(vzvarstr,&arg[9][2]);
    } else error->all("Illegal fix move command");

  } else error->all("Illegal fix move command");

  // optional args

  int scaleflag = 1;

  scaleflag=0;

  // error checks and warnings

  if (domain->dimension == 2) {
    if (mstyle == LINEAR && vzflag && vz != 0.0)
      error->all("Fix move cannot set linear z motion for 2d problem");
    if (mstyle == WIGGLE && azflag && az != 0.0)
      error->all("Fix move cannot set wiggle z motion for 2d problem");
    if (mstyle == RIGGLE && (axis[0] != 0.0 || axis[1] != 0.0))   
      error->all("Fix move cannot rotate aroung non z-axis for 2d problem");
    if (mstyle == ROTATE && (axis[0] != 0.0 || axis[1] != 0.0))
      error->all("Fix move cannot rotate aroung non z-axis for 2d problem");
    if (mstyle == VARIABLE && (zvarstr || vzvarstr))
      error->all("Fix move cannot define z or vz variable for 2d problem");
  }

  if (atom->angmom_flag && comm->me == 0)
    error->warning("Fix move does not update angular momentum");  
  if (atom->quat_flag && comm->me == 0)
    error->warning("Fix move does not update quaternions");       

  // setup scaling and apply scaling factors to velocity & amplitude

  if ((mstyle == LINEAR || mstyle == WIGGLE || mstyle == ROTATE || mstyle == RIGGLE || mstyle == DEFORM) && 
      scaleflag) {
    if (domain->lattice == NULL)
      error->all("Use of fix move with undefined lattice");

    double xscale,yscale,zscale;
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    if (mstyle == LINEAR) {
      if (vxflag) vx *= xscale;
      if (vyflag) vy *= yscale;
      if (vzflag) vz *= zscale;
    } else if (mstyle == WIGGLE) {
      if (axflag) ax *= xscale;
      if (ayflag) ay *= yscale;
      if (azflag) az *= zscale;
    } else if (mstyle == ROTATE || mstyle == RIGGLE || mstyle == DEFORM) {	
      point[0] *= xscale;
      point[1] *= yscale;
      point[2] *= zscale;
    }
  }

  // set omega_rotate from period

  if (mstyle == WIGGLE || mstyle == ROTATE || mstyle == RIGGLE || mstyle == DEFORM) { 
    double PI = 4.0 * atan(1.0);
    omega_rotate = 2.0*PI / period;
  }

  // runit = unit vector along rotation axis

  if (mstyle == ROTATE || mstyle == RIGGLE) { 
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all("Fix move cannot have 0 length rotation vector");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }

  // set omega_flag if particles store omega
  omega_flag = 0;  

  maxatom = 0; 

  displace = velocity = NULL;

  STL_tri= static_cast<FixMeshGran*>(modify->fix[f_i])->STLdata;
  STL_tri->initMove(skinSafetyFactor);

  x = STL_tri->x;  
  xoriginal=STL_tri->xoriginal; 
  v=STL_tri->v; 
  f=STL_tri->f; 
  nlocal = STL_tri->xvf_len; 
  rmass=STL_tri->rmass; 

  type=atom->type;
  mass=atom->mass;

  time_origin = update->ntimestep;

  restart_global = 1; 
}

/* ---------------------------------------------------------------------- */

FixMoveTri::~FixMoveTri()
{
  STL_tri->movingMesh=0;

  // delete locally stored arrays
  memory->destroy_2d_double_array(displace);
  memory->destroy_2d_double_array(velocity);

  delete [] xvarstr;
  delete [] yvarstr;
  delete [] zvarstr;
  delete [] vxvarstr;
  delete [] vyvarstr;
  delete [] vzvarstr;
}

/* ---------------------------------------------------------------------- */

int FixMoveTri::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveTri::init()
{
  dt = update->dt;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  // set indices and style of all variables

  if (mstyle == VARIABLE) {
    if (xvarstr) {
      xvar = input->variable->find(xvarstr);
      if (xvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(xvar)) xvarstyle = EQUAL;
      else if (input->variable->atomstyle(xvar)) xvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }
    if (yvarstr) {
      yvar = input->variable->find(yvarstr);
      if (yvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(yvar)) yvarstyle = EQUAL;
      else if (input->variable->atomstyle(yvar)) yvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }
    if (zvarstr) {
      zvar = input->variable->find(zvarstr);
      if (zvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(zvar)) zvarstyle = EQUAL;
      else if (input->variable->atomstyle(zvar)) zvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }
    if (vxvarstr) {
      vxvar = input->variable->find(vxvarstr);
      if (vxvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(vxvar)) vxvarstyle = EQUAL;
      else if (input->variable->atomstyle(vxvar)) vxvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }
    if (vyvarstr) {
      vyvar = input->variable->find(vyvarstr);
      if (vyvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(vyvar)) vyvarstyle = EQUAL;
      else if (input->variable->atomstyle(vyvar)) vyvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }
    if (vzvarstr) {
      vzvar = input->variable->find(vzvarstr);
      if (vzvar < 0) error->all("Variable name for fix move does not exist");
      if (input->variable->equalstyle(vzvar)) vzvarstyle = EQUAL;
      else if (input->variable->atomstyle(vzvar)) vzvarstyle = ATOM;
      else error->all("Variable for fix move is invalid style");
    }

    displaceflag = velocityflag = 0;
    if (xvarstr && xvarstyle == ATOM) displaceflag = 1;
    if (yvarstr && yvarstyle == ATOM) displaceflag = 1;
    if (zvarstr && zvarstyle == ATOM) displaceflag = 1;
    if (vxvarstr && vxvarstyle == ATOM) velocityflag = 1;
    if (vyvarstr && vyvarstyle == ATOM) velocityflag = 1;
    if (vzvarstr && vzvarstyle == ATOM) velocityflag = 1;
  }

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ----------------------------------------------------------------------
   set x,v of particles
------------------------------------------------------------------------- */

void FixMoveTri::initial_integrate(int vflag)
{
  double dtfm;
  double xold[3],a[3],b[3],c[3],d[3],disp[3];
  double ddotr,dx,dy,dz;

  double delta = (update->ntimestep - time_origin) * dt;

  STL_tri->vecToPoint();

  // for linear: X = X0 + V*dt

  if (mstyle == LINEAR) {
    for (int i = 0; i < nlocal; i++) {
    
	xold[0] = x[i][0];
	xold[1] = x[i][1];
	xold[2] = x[i][2];

	if (vxflag) {
	  v[i][0] = vx;
	  x[i][0] = xoriginal[i][0] + vx*delta;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  x[i][0] += dtv * v[i][0];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][0] += dtfm * f[i][0];
	  x[i][0] += dtv * v[i][0];
	}

	if (vyflag) {
	  v[i][1] = vy;
	  x[i][1] = xoriginal[i][1] + vy*delta;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][1] += dtfm * f[i][1];
	  x[i][1] += dtv * v[i][1];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][1] += dtfm * f[i][1];
	  x[i][1] += dtv * v[i][1];
	}

	if (vzflag) {
	  v[i][2] = vz;
	  x[i][2] = xoriginal[i][2] + vz*delta;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][2] += dtfm * f[i][2];
	  x[i][2] += dtv * v[i][2];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][2] += dtfm * f[i][2];
	  x[i][2] += dtv * v[i][2];
	}

    }

  // for wiggle: X = X0 + A sin(w*dt)

  } else if (mstyle == WIGGLE) {
    double arg = omega_rotate * delta;
    double sine = sin(arg);
    double cosine = cos(arg);

    for (int i = 0; i < nlocal; i++) {
      
	xold[0] = x[i][0];
	xold[1] = x[i][1];
	xold[2] = x[i][2];

	if (axflag) {
	  v[i][0] = ax*omega_rotate*cosine;
	  x[i][0] = xoriginal[i][0] + ax*sine;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  x[i][0] += dtv * v[i][0];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][0] += dtfm * f[i][0];
	  x[i][0] += dtv * v[i][0];
	}

	if (ayflag) {
	  v[i][1] = ay*omega_rotate*cosine;
	  x[i][1] = xoriginal[i][1] + ay*sine;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][1] += dtfm * f[i][1];
	  x[i][1] += dtv * v[i][1];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][1] += dtfm * f[i][1];
	  x[i][1] += dtv * v[i][1];
	}

	if (azflag) {
	  v[i][2] = az*omega_rotate*cosine;
	  x[i][2] = xoriginal[i][2] + az*sine;
	} else if (rmass) {
	  dtfm = dtf / rmass[i];
	  v[i][2] += dtfm * f[i][2];
	  x[i][2] += dtv * v[i][2];
	} else {
	  dtfm = dtf / mass[type[i]];
	  v[i][2] += dtfm * f[i][2];
	  x[i][2] += dtv * v[i][2];
	}

    }

  // for rotate by right-hand rule around omega:
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // w = omega of rotation (from period)
  // X0 = xoriginal = initial coord of atom
  // R0 = runit = unit vector for R
  // D = X0 - P = vector from P to X0
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X0
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(w*dt) + B sin(w*dt)
  // V = w R0 cross (A cos(w*dt) + B sin(w*dt))

  } else if (mstyle == ROTATE || mstyle == RIGGLE) {
    
    double arg;

    if (mstyle == ROTATE) arg = omega_rotate * delta;
    if (mstyle == RIGGLE)
    {
        arg = ampl*sin(omega_rotate * delta);
        
    }
    
    double sine = sin(arg);
    double cosine = cos(arg);

    for (int i = 0; i < nlocal; i++) {
      
	xold[0] = x[i][0];
	xold[1] = x[i][1];
	xold[2] = x[i][2];

	d[0] = xoriginal[i][0] - point[0];
	d[1] = xoriginal[i][1] - point[1];
	d[2] = xoriginal[i][2] - point[2];
	ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
	c[0] = ddotr*runit[0];
	c[1] = ddotr*runit[1];
	c[2] = ddotr*runit[2];
	a[0] = d[0] - c[0];
	a[1] = d[1] - c[1];
	a[2] = d[2] - c[2];
	b[0] = runit[1]*a[2] - runit[2]*a[1];
	b[1] = runit[2]*a[0] - runit[0]*a[2];
	b[2] = runit[0]*a[1] - runit[1]*a[0];
	disp[0] = a[0]*cosine  + b[0]*sine;
	disp[1] = a[1]*cosine  + b[1]*sine;
	disp[2] = a[2]*cosine  + b[2]*sine;

	x[i][0] = point[0] + c[0] + disp[0];
	x[i][1] = point[1] + c[1] + disp[1];
	x[i][2] = point[2] + c[2] + disp[2];

	if (mstyle == ROTATE)
	{
        v[i][0] = omega_rotate * (runit[1]*disp[2] - runit[2]*disp[1]);
        v[i][1] = omega_rotate * (runit[2]*disp[0] - runit[0]*disp[2]);
        v[i][2] = omega_rotate * (runit[0]*disp[1] - runit[1]*disp[0]);
	}
	if (mstyle == RIGGLE)
	{
        v[i][0] = ampl*cos(omega_rotate * delta)/omega_rotate * (runit[1]*disp[2] - runit[2]*disp[1]);
        v[i][1] = ampl*cos(omega_rotate * delta)/omega_rotate * (runit[2]*disp[0] - runit[0]*disp[2]);
        v[i][2] = ampl*cos(omega_rotate * delta)/omega_rotate * (runit[0]*disp[1] - runit[1]*disp[0]);
	}

    }

  } else if (mstyle == DEFORM ) { 
    error->all("need to account for length change in ogK");

    double scaling = ampl*sin(omega_rotate * delta);
    double time_deriv = ampl*cos(omega_rotate * delta)/omega_rotate;

    for (int i = 0; i < nlocal; i++) {
        x[i][0] = xoriginal[i][0] + scaling*(xoriginal[i][0]-point[0]);
        x[i][1] = xoriginal[i][1] + scaling*(xoriginal[i][1]-point[1]);
        x[i][2] = xoriginal[i][2] + scaling*(xoriginal[i][2]-point[2]);

        v[i][0] = time_deriv * (xoriginal[i][0]-point[0]);
        v[i][1] = time_deriv * (xoriginal[i][1]-point[1]);
        v[i][2] = time_deriv * (xoriginal[i][2]-point[2]);
    }

  // for variable: compute x,v from variables
  } else if (mstyle == VARIABLE) {

    // reallocate displace and velocity arrays as necessary

    if ((displaceflag || velocityflag) && nlocal > maxatom) {
      maxatom = nlocal; 
      if (displaceflag) {
	memory->destroy_2d_double_array(displace);
	displace = memory->create_2d_double_array(maxatom,3,"move:displace");
      }
      if (velocityflag) {
	memory->destroy_2d_double_array(velocity);
	velocity = memory->create_2d_double_array(maxatom,3,"move:velocity");
      }
    }

    // pre-compute variable values

    if (xvarstr) {
      if (xvarstyle == EQUAL) dx = input->variable->compute_equal(xvar);
      else input->variable->compute_atom(xvar,igroup,&displace[0][0],3,0);
    }
    if (yvarstr) {
      if (yvarstyle == EQUAL) dy = input->variable->compute_equal(yvar);
      else input->variable->compute_atom(yvar,igroup,&displace[0][1],3,0);
    }
    if (zvarstr) {
      if (zvarstyle == EQUAL) dz = input->variable->compute_equal(zvar);
      else input->variable->compute_atom(zvar,igroup,&displace[0][2],3,0);
    }
    if (vxvarstr) {
      if (vxvarstyle == EQUAL) vx = input->variable->compute_equal(vxvar);
      else input->variable->compute_atom(vxvar,igroup,&velocity[0][0],3,0);
    }
    if (vyvarstr) {
      if (vyvarstyle == EQUAL) vy = input->variable->compute_equal(vyvar);
      else input->variable->compute_atom(vyvar,igroup,&velocity[0][1],3,0);
    }
    if (vzvarstr) {
      if (vzvarstyle == EQUAL) vz = input->variable->compute_equal(vzvar);
      else input->variable->compute_atom(vzvar,igroup,&velocity[0][2],3,0);
    }

    // update x,v

    for (int i = 0; i < nlocal; i++) {
      
	xold[0] = x[i][0];
	xold[1] = x[i][1];
	xold[2] = x[i][2];

	if (xvarstr && vxvarstr) {
	  if (vxvarstyle == EQUAL) v[i][0] = vx;
	  else v[i][0] = velocity[i][0];
	  if (xvarstyle == EQUAL) x[i][0] = xoriginal[i][0] + dx;
	  else x[i][0] = xoriginal[i][0] + displace[i][0];
	} else if (xvarstr) {
	  if (xvarstyle == EQUAL) x[i][0] = xoriginal[i][0] + dx;
	  else x[i][0] = xoriginal[i][0] + displace[i][0];
	} else if (vxvarstr) {
	  if (vxvarstyle == EQUAL) v[i][0] = vx;
	  else v[i][0] = velocity[i][0];
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    x[i][0] += dtv * v[i][0];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    x[i][0] += dtv * v[i][0];
	  }
	} else {
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    v[i][0] += dtfm * f[i][0];
	    x[i][0] += dtv * v[i][0];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    v[i][0] += dtfm * f[i][0];
	    x[i][0] += dtv * v[i][0];
	  }
	}

	if (yvarstr && vyvarstr) {
	  if (vyvarstyle == EQUAL) v[i][1] = vy;
	  else v[i][1] = velocity[i][1];
	  if (yvarstyle == EQUAL) x[i][1] = xoriginal[i][1] + dy;
	  else x[i][1] = xoriginal[i][1] + displace[i][1];
	} else if (yvarstr) {
	  if (yvarstyle == EQUAL) x[i][1] = xoriginal[i][1] + dy;
	  else x[i][1] = xoriginal[i][1] + displace[i][1];
	} else if (vyvarstr) {
	  if (vyvarstyle == EQUAL) v[i][1] = vy;
	  else v[i][1] = velocity[i][1];
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    x[i][1] += dtv * v[i][1];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    x[i][1] += dtv * v[i][1];
	  }
	} else {
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    v[i][1] += dtfm * f[i][1];
	    x[i][1] += dtv * v[i][1];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    v[i][1] += dtfm * f[i][1];
	    x[i][1] += dtv * v[i][1];
	  }
	}

	if (zvarstr && vzvarstr) {
	  if (vzvarstyle == EQUAL) v[i][2] = vz;
	  else v[i][2] = velocity[i][2];
	  if (zvarstyle == EQUAL) x[i][2] = xoriginal[i][2] + dz;
	  else x[i][2] = xoriginal[i][2] + displace[i][2];
	} else if (zvarstr) {
	  if (zvarstyle == EQUAL) x[i][2] = xoriginal[i][2] + dz;
	  else x[i][2] = xoriginal[i][2] + displace[i][2];
	} else if (vzvarstr) {
	  if (vzvarstyle == EQUAL) v[i][2] = vz;
	  else v[i][2] = velocity[i][2];
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    x[i][2] += dtv * v[i][2];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    x[i][2] += dtv * v[i][2];
	  }
	} else {
	  if (rmass) {
	    dtfm = dtf / rmass[i];
	    v[i][2] += dtfm * f[i][2];
	    x[i][2] += dtv * v[i][2];
	  } else {
	    dtfm = dtf / mass[type[i]];
	    v[i][2] += dtfm * f[i][2];
	    x[i][2] += dtv * v[i][2];
	  }
	}

    }
  }

  //endpoints are transformed back into vectors 
  STL_tri->pointToVec();
}

/* ----------------------------------------------------------------------
   final NVE of particles with NULL components
------------------------------------------------------------------------- */

void FixMoveTri::final_integrate()
{
    
}

/* ---------------------------------------------------------------------- */

void FixMoveTri::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  // outermost level - update v and x
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMoveTri::final_integrate_respa(int ilevel)
{
  if (ilevel == nlevels_respa-1) final_integrate();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMoveTri::write_restart(FILE *fp)
{
  int n = 0;
  int listlen=1+STL_tri->xvf_len*3;
  double *list=new double[listlen];
  list[n++] = time_origin;

  STL_tri->vecToPoint();
  for (int i=0;i<STL_tri->xvf_len;i++)
  {
     for(int j=0;j<3;j++) list[n++]=STL_tri->xoriginal[i][j];
  }
  STL_tri->pointToVec();

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
  delete []list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMoveTri::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  time_origin = static_cast<int> (list[n++]);

  for (int i=0;i<STL_tri->xvf_len;i++)
  {
       for(int j=0;j<3;j++) STL_tri->xoriginal[i][j]=list[n++];
  }
}
