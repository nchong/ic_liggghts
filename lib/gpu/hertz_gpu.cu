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

#include "cuPrintf.cu"
#include "fshearmap.cu"

static NVCTimer timer;

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

#ifdef VERBOSE
  printf("> hertz_gpu_init\n");
  printf("> cell_size = %f\n", cell_size);
  printf("> boxhi = {%f, %f, %f}; boxlo = {%f, %f, %f}\n", 
    boxhi[0], boxhi[1], boxhi[2],
    boxlo[0], boxlo[1], boxlo[2]);
  printf("> ncellx = %d; ncelly = %d; ncellz = %d;\n",
    ncellx, ncelly, ncellz);
#endif
   
  init_cell_list_const(cell_size, skin, boxlo, boxhi);

  timer.init();

  return true;
}

// ---------------------------------------------------------------------------
// Clear memory on host and device
// ---------------------------------------------------------------------------
EXTERN void hertz_gpu_clear() {
  clear_cell_list(cell_list_gpu);
  ASSERT_NO_CUDA_ERROR(cudaFree(d_x));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_v));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_omega));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_torque));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_f));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_atom_type));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_radius));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_rmass));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_Yeff));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_Geff));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_betaeff));
  ASSERT_NO_CUDA_ERROR(cudaFree(d_coeffFrict));
}

EXTERN struct fshearmap *create_fshearmap(
  //number of atoms (ranged over by ii)
  int inum, /*list->inum*/
  //id of atom ii (ranged over by i)
  int *ilist, /*list->ilist*/
  //numneigh[i] is the number of neighbor atoms to atom i (ranged over by jj)
  int *numneigh, /*list->numneigh*/
  //firstneigh[i][jj] is the id of atom jj (ranged over by j)
  int **firstneigh, /*list->firstneigh*/
  //firsttouch[i][jj] = 1, if i and j are in contact
  //                    0, otherwise
  int **firsttouch, /*listgranhistory->firstneigh*/
  //firstshear[i][3*jj] is the shear vector between atoms i and j
  double **firstshear /*listgranhistory->firstdouble*/) {
  struct fshearmap *map = malloc_fshearmap(inum);

  for (int ii=0; ii<inum; ii++) {
    int i = ilist[ii];
    assert(i < inum);
    int jnum = numneigh[i];
    assert(jnum < 32);

    for (int jj = 0; jj<jnum; jj++) {
      int j = firstneigh[i][jj];

      //TODO: necessary to check firsttouch[i][jj] == 1?
      double *shear = &firstshear[i][3*jj];
      insert_fshearmap(map, i, j, shear);

      //symmetric-shear for particle j
      double nshear[3];
      nshear[0] = -shear[0];
      nshear[1] = -shear[1];
      nshear[2] = -shear[2];
      insert_fshearmap(map, j, i, nshear);
    }
  }

  //paranoid
  for (int ii=0; ii<inum; ii++) {
    int i = ilist[ii];
    int jnum = numneigh[i];

    for (int jj = 0; jj<jnum; jj++) {
      int j = firstneigh[i][jj];

      double *shear = &firstshear[i][3*jj];
      double *result = retrieve_fshearmap(map, i, j);
      assert(result);
      assert(result[0] == shear[0]);
      assert(result[1] == shear[1]);
      assert(result[2] == shear[2]);

      result = retrieve_fshearmap(map, j, i);
      assert(result);
      assert(result[0] == -shear[0]);
      assert(result[1] == -shear[1]);
      assert(result[2] == -shear[2]);
    }
  }

  return map;
}

EXTERN void update_from_fshearmap(
  struct fshearmap *map,
  int inum,
  int *ilist,
  int *numneigh,
  int **firstneigh,
  int **firsttouch,
  double **firstshear) {
  for (int ii=0; ii<inum; ii++) {
    int i = ilist[ii];
    int jnum = numneigh[i];

    for (int jj = 0; jj<jnum; jj++) {
      int j = firstneigh[i][jj];

      double *result = retrieve_fshearmap(map, i, j);
      assert(result);
      double *shear = &firstshear[i][3*jj];
      shear[0] = result[0];
      shear[1] = result[1];
      shear[2] = result[2];
      //TODO: necessary to uncheck firsttouch[i][jj]?
    }
  }
  free_fshearmap(map);
}

EXTERN double hertz_gpu_cell(
  const bool eflag, const bool vflag,
  const int inum, const int nall, const int ago,

  //inputs
  double **host_x, double **host_v, double **host_omega,
  double *host_radius, double *host_rmass,
  int *host_type,
  double dt,

  //stiffness parameters
  int num_atom_types,
  double **host_Yeff,
  double **host_Geff,
  double **host_betaeff,
  double **host_coeffFrict,
  double nktv2p,

  //inouts 
  struct fshearmap *&host_fshearmap,
  double **host_torque, double **host_force)
{
  static int blockSize = BLOCK_1D;
  static int ncell = ncellx*ncelly*ncellz;

#ifdef VERBOSE
  printf("> inum = %d\n", inum);
  printf("> nall = %d\n", nall);
#endif
  // -----------------------------------------------------------
  // Device variables
  // -----------------------------------------------------------
  const int SIZE_1D = (nall * sizeof(double));
  const int SIZE_2D = (3 * SIZE_1D);
  static bool first_call = true;
  if (first_call) {
    first_call = false;
    printf("> initialising device datastructures\n");

    //malloc device data
    init_cell_list(cell_list_gpu, nall, ncell, blockSize);
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_x, SIZE_2D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_v, SIZE_2D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_omega, SIZE_2D));
    //shear done by malloc_device_fshearmap()
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_torque, SIZE_2D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_f, SIZE_2D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_atom_type, SIZE_1D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_radius, SIZE_1D));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_rmass, SIZE_1D));

    const int PARAM_SIZE = num_atom_types * num_atom_types * sizeof(double);
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_Yeff, PARAM_SIZE));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_Geff, PARAM_SIZE));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_betaeff, PARAM_SIZE));
    ASSERT_NO_CUDA_ERROR(cudaMalloc((void**)&d_coeffFrict, PARAM_SIZE));

    //we flatten atom vectors via these host arrays
    h_x = (double *)malloc(SIZE_2D);
    h_v = (double *)malloc(SIZE_2D);
    h_omega = (double *)malloc(SIZE_2D);
    h_torque = (double *)malloc(SIZE_2D);
    h_f = (double *)malloc(SIZE_2D);

    assert(h_x);
    assert(h_v);
    assert(h_omega);
    assert(h_torque);
    assert(h_f);

    //flatten stiffness lookup tables
    double *h_Yeff = (double *)malloc(PARAM_SIZE);
    double *h_Geff = (double *)malloc(PARAM_SIZE);
    double *h_betaeff = (double *)malloc(PARAM_SIZE);
    double *h_coeffFrict = (double *)malloc(PARAM_SIZE);

    assert(h_Yeff);
    assert(h_Geff);
    assert(h_betaeff);
    assert(h_coeffFrict);

    for (int i=0; i<num_atom_types; i++) {
      for (int j=0; j<num_atom_types; j++) {
        h_Yeff[i + (j*num_atom_types)] = host_Yeff[i][j];
        h_Geff[i + (j*num_atom_types)] = host_Geff[i][j];
        h_betaeff[i + (j*num_atom_types)] = host_betaeff[i][j];
        h_coeffFrict[i + (j*num_atom_types)] = host_coeffFrict[i][j];
      }
    }

    //stiffness coefficients initialised once
    ASSERT_NO_CUDA_ERROR(
      cudaMemcpy(d_Yeff, h_Yeff, PARAM_SIZE, cudaMemcpyHostToDevice));
    ASSERT_NO_CUDA_ERROR(
      cudaMemcpy(d_Geff, h_Geff, PARAM_SIZE, cudaMemcpyHostToDevice));
    ASSERT_NO_CUDA_ERROR(
      cudaMemcpy(d_betaeff, h_betaeff, PARAM_SIZE, cudaMemcpyHostToDevice));
    ASSERT_NO_CUDA_ERROR(
      cudaMemcpy(d_coeffFrict, h_coeffFrict, PARAM_SIZE, cudaMemcpyHostToDevice));

    free(h_Yeff);
    free(h_Geff);
    free(h_betaeff);
    free(h_coeffFrict);
  }

  //flatten incoming 2d atom data
  for (int i=0; i< nall; i++) {
    for (int j=0; j<3; j++) {
      h_x[(i*3)+j] = host_x[i][j];
      h_v[(i*3)+j] = host_v[i][j];
      h_omega[(i*3)+j] = host_omega[i][j];
      h_torque[(i*3)+j] = host_torque[i][j];
      h_f[(i*3)+j] = host_force[i][j];
    }
  }

  //copy across 2d atom data
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_x, h_x, SIZE_2D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_v, h_v, SIZE_2D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_omega, h_omega, SIZE_2D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_torque, h_torque, SIZE_2D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_f, h_f, SIZE_2D, cudaMemcpyHostToDevice));

  //just copy across 1d atom data
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_atom_type, host_type, SIZE_1D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_radius, host_radius, SIZE_1D, cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_rmass, host_rmass, SIZE_1D, cudaMemcpyHostToDevice));

  struct fshearmap gpu_fshearmap;
  malloc_device_fshearmap(gpu_fshearmap, inum);
  copy_into_device_fshearmap(gpu_fshearmap, host_fshearmap);

  build_cell_list(host_x[0], host_type, cell_list_gpu, 
	  ncell, ncellx, ncelly, ncellz, blockSize, inum, nall, ago);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Pre-kernel error: %s.\n", cudaGetErrorString(err));
    exit(1);
  }

  cudaPrintfInit();

  timer.start();
  const int BX=blockSize;
  dim3 GX(ncellx, ncelly*ncellz);
  kernel_hertz_cell<true,true,64><<<GX,BX>>>(
    cell_list_gpu.pos,
    cell_list_gpu.idx,
    cell_list_gpu.type,
    cell_list_gpu.natom,
    inum, ncellx, ncelly, ncellz,
    d_x, d_v, d_omega, 
    d_atom_type, d_radius, d_rmass, 

    gpu_fshearmap.valid, gpu_fshearmap.key, gpu_fshearmap.shear,
    d_torque, d_f,
    dt, num_atom_types,
    d_Yeff, d_Geff, d_betaeff, d_coeffFrict, nktv2p
  );
  timer.stop();
  timer.add_to_total();

  cudaThreadSynchronize();
  cudaPrintfDisplay(stderr, true);
  cudaPrintfEnd();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("Post-kernel error: %s.\n", cudaGetErrorString(err));
    exit(1);
  }

  //copy back force calculations (shear, torque, force)
  copy_from_device_fshearmap(gpu_fshearmap, host_fshearmap);
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(h_torque, d_torque, SIZE_2D, cudaMemcpyDeviceToHost));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(h_f, d_f, SIZE_2D, cudaMemcpyDeviceToHost));
  for (int i=0; i<nall; i++) {
    host_torque[i][0] = h_torque[(i*3)];
    host_torque[i][1] = h_torque[(i*3)+1];
    host_torque[i][2] = h_torque[(i*3)+2];

    host_force[i][0] = h_f[(i*3)];
    host_force[i][1] = h_f[(i*3)+1];
    host_force[i][2] = h_f[(i*3)+2];
  }

  //now free all malloc'ed data
  free_device_fshearmap(gpu_fshearmap);
  //gpu shearmap freed by call to update_from_fshearmap
  //all other atom data is static

  return 0.0;
}

EXTERN void hertz_gpu_time() {
  printf("Kernel time %.16fms\n", timer.total_time());
}

EXTERN int hertz_gpu_num_devices() {
  return 0;
}

EXTERN double hertz_gpu_bytes() {
  return 0.0;
}
