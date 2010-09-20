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

#ifndef LMP_PARTICLE_TO_INSERT_H
#define LMP_PARTICLE_TO_INSERT_H

#include "memory.h"
#include "pointers.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsert : protected Pointers
    {
      public:
        ParticleToInsert(LAMMPS* lmp,int ns);
        ~ParticleToInsert();

        void random_rotate(double,double,double);

        int nspheres;
        int atom_type;
        double density_ins;
        double volume_ins;
        double mass_ins;
        double *radius_ins;
        double r_bound;
        double **x_ins;

        double *xcm;
        double *inertia;
        double *ex_space,*ey_space,*ez_space;
        double **displace;

    };

}

#endif

