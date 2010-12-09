#ifndef CUDA_FSHEARMAP
#define CUDA_FSHEARMAP

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.cuh"

#define ILEN 32 //<number of elements per i slot

//---------------------------------------------------------------------------
// Datastructure
//---------------------------------------------------------------------------

struct fshearmap {
  int inum; //number of i slots
  bool *valid;
  int *key;
  double *shear;
};

inline int hash_of(int i, int key) {
  int hash = (i*ILEN) + (key%ILEN);
  const int lower = (i * ILEN) - 1;
  const int upper = ((i+1) * ILEN);
  assert(lower < hash && hash < upper);
  return hash;
}

//---------------------------------------------------------------------------
// Host functions
//---------------------------------------------------------------------------
struct fshearmap *malloc_fshearmap(int inum) {
  struct fshearmap *map = 
    (struct fshearmap *)malloc(sizeof(struct fshearmap));
  assert(map);
  map->inum = inum;

  const int TABLE_SIZE = inum * ILEN;
  map->valid = (bool *)malloc(TABLE_SIZE * sizeof(bool));
  assert(map->valid);
  memset(map->valid, false, TABLE_SIZE * sizeof(bool));

  map->key = (int *)malloc(TABLE_SIZE * sizeof(int));
  assert(map->key);
  memset(map->key, 0, TABLE_SIZE * sizeof(int));

  map->shear = (double *)malloc(TABLE_SIZE * 3 * sizeof(double));
  assert(map->shear);
  memset(map->shear, 0, TABLE_SIZE * 3 * sizeof(double));

  return map;
}

void free_fshearmap(struct fshearmap *map) {
  free(map->valid);  
  free(map->key);  
  free(map->shear);  
  free(map);
}

double *retrieve_fshearmap(struct fshearmap *map, int i, int key) {
  double *result = NULL;

  int hash = hash_of(i, key);
  if (map->valid[hash] && map->key[hash] == key) {
    return &map->shear[hash*3];
  } else {
    int probe = hash;
    do {
      probe = hash_of(i, probe+1);
    } while (map->key[probe] != key && probe != hash);

    if (map->valid[probe] && map->key[probe] == key && probe != hash) {
      return &map->shear[probe*3];
    }
  }
  return result;
}

void insert_fshearmap(struct fshearmap *&map, int i, int key, double shear[3]) {
  assert(0<=i && i<map->inum);

  //linear probe on collision
  int hash = hash_of(i, key);
  if (map->valid[hash]) {
    int probe = hash;
    do {
      probe = hash_of(i, probe+1);
    } while (map->valid[probe] && probe != hash);

    if (probe == hash) {
      printf("error: map is full!\n"); //could resize here
      exit(1);
    } else {
      hash = probe;
    }
  }

  assert(map->valid[hash] == false);
  map->valid[hash] = true;
  map->key[hash] = key;
  map->shear[(hash*3)] = shear[0];
  map->shear[(hash*3)+1] = shear[1];
  map->shear[(hash*3)+2] = shear[2];
}

void print_fshearmap(struct fshearmap *map) {
  for (int i=0; i<map->inum; i++) {
    printf("map[%i]\n", i);
    for (int j=0; j<ILEN; j++) {
      int hash = hash_of(i, j);
      printf("%d: ", hash);
      if (map->valid[hash]) {
        printf("%d |-> {%f, %f, %f}\n",
          map->key[hash],
          map->shear[(hash*3)],
          map->shear[(hash*3)+1],
          map->shear[(hash*3)+2]);
      } else {
        printf("#\n");
      }
    }
  }
}

//---------------------------------------------------------------------------
// Device functions
//---------------------------------------------------------------------------

// NB: gpu_map is a *host* side datastructure that contains device pointers
void malloc_device_fshearmap(struct fshearmap &gpu_map, int inum) {
  const int TABLE_SIZE = inum * ILEN;
  ASSERT_NO_CUDA_ERROR(
    cudaMalloc((void **)&(gpu_map.valid), TABLE_SIZE * sizeof(bool)));
  ASSERT_NO_CUDA_ERROR(
    cudaMalloc((void **)&(gpu_map.key), TABLE_SIZE * sizeof(int)));
  ASSERT_NO_CUDA_ERROR(
    cudaMalloc((void **)&(gpu_map.shear), TABLE_SIZE * 3 * sizeof(double)));
}

void free_device_fshearmap(struct fshearmap &gpu_map) {
  ASSERT_NO_CUDA_ERROR(cudaFree(gpu_map.valid));
  ASSERT_NO_CUDA_ERROR(cudaFree(gpu_map.key));
  ASSERT_NO_CUDA_ERROR(cudaFree(gpu_map.shear));
}

void copy_into_device_fshearmap(struct fshearmap &gpu_map, struct fshearmap *map) {
  const int TABLE_SIZE = map->inum * ILEN;
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(gpu_map.valid, map->valid, TABLE_SIZE * sizeof(bool), cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(gpu_map.key, map->key, TABLE_SIZE * sizeof(int), cudaMemcpyHostToDevice));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(gpu_map.shear, map->shear, TABLE_SIZE * 3 * sizeof(double), cudaMemcpyHostToDevice));
}

void copy_from_device_fshearmap(struct fshearmap &gpu_map, struct fshearmap *map) {
  const int TABLE_SIZE = map->inum * ILEN;
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(map->valid, gpu_map.valid, TABLE_SIZE * sizeof(bool), cudaMemcpyDeviceToHost));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(map->key, gpu_map.key, TABLE_SIZE * sizeof(int), cudaMemcpyDeviceToHost));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(map->shear, gpu_map.shear, TABLE_SIZE * 3 * sizeof(double), cudaMemcpyDeviceToHost));
}

#define CUDA_HASH_OF(i, key) ((i*ILEN) + (key%ILEN))
__device__ double *cuda_retrieve_fshearmap(
  bool *map_valid, int *map_key, double *map_shear,
  int i, int key) {
  double *result = NULL;

  int hash = CUDA_HASH_OF(i, key);
  if (map_valid[hash] && map_key[hash] == key) {
    return &map_shear[hash*3];
  } else {
    int probe = hash;
    do {
      probe = CUDA_HASH_OF(i, (probe+1));
    } while (map_key[probe] != key && probe != hash);

    if (map_valid[probe] && map_key[probe] == key && probe != hash) {
      return &map_shear[probe*3];
    }
  }
  return result;
}

#endif
