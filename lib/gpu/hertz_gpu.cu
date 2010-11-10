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

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#include <iostream>
#include <cassert>
#include "nvc_macros.h"
#include "nvc_timer.h"
#include "nvc_device.h"
#include "pair_gpu_texture.h"
#include "pair_gpu_cell.h"
#include "lj_gpu_memory.cu"
#include "hertz_gpu_kernel.h"

// ---------------------------------------------------------------------------
// Return string with GPU info
// ---------------------------------------------------------------------------
EXTERN void hertz_gpu_name(const int id, const int max_nbors, char * name) {
  strcpy(name,"hertz");
}

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
EXTERN bool hertz_gpu_init(
  double *boxlo, double *boxhi, 
	double cell_size, double skin)
{

  ncellx = ceil(((boxhi[0] - boxlo[0]) + 2.0*cell_size) / cell_size);
  ncelly = ceil(((boxhi[1] - boxlo[1]) + 2.0*cell_size) / cell_size);
  ncellz = ceil(((boxhi[2] - boxlo[2]) + 2.0*cell_size) / cell_size);

  printf("> hertz_gpu_init\n");
  printf("> cell_size = %f\n", cell_size);
  printf("> boxhi = {%f, %f, %f}; boxlo = {%f, %f, %f}\n", 
    boxhi[0], boxhi[1], boxhi[2],
    boxlo[0], boxlo[1], boxlo[2]);
  printf("> ncellx = %d; ncelly = %d; ncellz = %d;\n",
    ncellx, ncelly, ncellz);
   
  init_cell_list_const(cell_size, skin, boxlo, boxhi);

  return true;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
EXTERN void hertz_gpu_clear() {
}

template <class numtyp, class acctyp>
double _hertz_gpu_cell()
{
  return 0.0;
}

EXTERN double hertz_gpu_cell()
{
  return 0.0;
}

EXTERN void hertz_gpu_time() {
}

EXTERN int hertz_gpu_num_devices() {
  return 0;
}

EXTERN double hertz_gpu_bytes() {
  return 0.0;
}
