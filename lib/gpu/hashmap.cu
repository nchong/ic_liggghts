
#ifndef CUDA_HASHMAP
#define CUDA_HASHMAP

#define ASSERT_NO_CUDA_ERROR( callReturningErrorstatus ) {     \
  cudaError_t err = callReturningErrorstatus;                  \
  if (err != cudaSuccess) {                                    \
    fprintf(stderr,                                            \
            "Cuda error (%s/%d) in file '%s' in line %i\n",    \
            cudaGetErrorString(err), err, __FILE__, __LINE__); \
    exit(1);                                                   \
  }                                                            \
} while(0);

#include <stdio.h>
#include <assert.h>

// --------------------------------------------------------------------------
// Hashmap datastructure
// --------------------------------------------------------------------------
struct entry {
  bool valid;
  int  key;
  double shear[3]; //data
};

struct hashmap {
  int size;
  struct entry *map;
};

// --------------------------------------------------------------------------
// C helper functions
// --------------------------------------------------------------------------
struct hashmap *create_hashmap(int size) {
  struct hashmap *hm = (struct hashmap *) malloc(sizeof(struct hashmap));
  if (!hm) {
    printf("error: could not malloc hashmap struct\n");
    exit(1);
  }
  hm->size = size;
  hm->map = (struct entry *) malloc(sizeof(struct entry) * size);
  if (!hm->map) {
    printf("error: could not malloc map of size %d\n", size);
    exit(1);
  }
  for (int i=0; i<size; i++) {
    hm->map[i].valid = false;
    hm->map[i].key = 0;
    hm->map[i].shear[0] = 0;
    hm->map[i].shear[1] = 0;
    hm->map[i].shear[2] = 0;
  }
  assert(hm);
  assert(hm->size == size);
  assert(hm->map);
  return hm;
}

void free_hashmap(struct hashmap *hm) {
  free(hm->map);
  free(hm);
}

void free_device_hashmap(struct hashmap *d_hm) {
  struct entry *d_map;
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(&d_map, &d_hm->map, sizeof(&d_map), cudaMemcpyDeviceToHost));
  ASSERT_NO_CUDA_ERROR(
    cudaFree(d_map)
  );
  ASSERT_NO_CUDA_ERROR(
    cudaFree(d_hm)
  );
}

struct entry *insert_hashmap(struct hashmap *hm, int key, double shear[3]) {
  assert(hm->size > 0);
  int hash = key % (hm->size);

  //linear probe on collision
  if (hm->map[hash].valid) {
    int probe = (hash + 1) % hm->size;
    while (hm->map[probe].valid && probe != hash) {
      probe = (probe + 1) % hm->size;
    }
    
    if (probe == hash) {
      printf("error: hashmap is full!\n"); //could resize here
      exit(1);
    } else {
      hash = probe;
    }
  }

  //straightforward insertion
  hm->map[hash].valid    = true;
  hm->map[hash].key      = key;
  hm->map[hash].shear[0] = shear[0];
  hm->map[hash].shear[1] = shear[1];
  hm->map[hash].shear[2] = shear[2];
  return &hm->map[hash];
}

struct entry *retrieve_hashmap(struct hashmap *hm, int key) {
  assert(hm);
  assert(hm->map);
  assert(hm->size > 0 );

  struct entry *result = NULL;
  int hash = key % (hm->size);
  if (hm->map[hash].valid && hm->map[hash].key == key) {
    result = &hm->map[hash];
  } else {
    int probe = (hash + 1) % hm->size;
    while (hm->map[probe].key != key && probe != hash) {
      probe = (probe + 1) % hm->size;
    }

    if (probe != hash && hm->map[probe].valid) {
      result = &hm->map[probe];
    }
  }
  return result;
}

void print_hashmap(struct hashmap *hm) {
  for (int i=0; i<hm->size; i++) {
    printf("%d: ", i);
    if (hm->map[i].valid) {
      printf("%d |-> {%f, %f, %f}\n", 
        hm->map[i].key,
        hm->map[i].shear[0],
        hm->map[i].shear[1],
        hm->map[i].shear[2]);
    } else {
      printf("#\n");
    }
  }
}

/*
 * Create a device copy of a C hashmap
 */
struct hashmap *c_to_device_hashmap(struct hashmap *hm) {
  struct entry *d_map;
  size_t table_size = sizeof(struct entry) * hm->size;
  ASSERT_NO_CUDA_ERROR(
    cudaMalloc((void **)&d_map, table_size));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(d_map, hm->map, table_size, cudaMemcpyHostToDevice));

  struct hashmap *d_hm;
  ASSERT_NO_CUDA_ERROR(
    cudaMalloc((void **)&d_hm, sizeof(struct hashmap)));

  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(&d_hm->size, &hm->size, sizeof(int), cudaMemcpyHostToDevice));
  //nb copying sizeof pointer *not* the table here
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(&d_hm->map, &d_map, sizeof(&d_map), cudaMemcpyHostToDevice));

  return d_hm;
}

/*
 * Inverse of c_to_device_hashmap. Create a C copy of a device hashmap
 */
struct hashmap *device_to_c_hashmap(struct hashmap *d_hm) {
  struct hashmap *hm = (struct hashmap *) malloc(sizeof(struct hashmap));
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(&hm->size, &d_hm->size, sizeof(int), cudaMemcpyDeviceToHost));

  //nb this is a *device* pointer to unwrap the hashmap struct
  //doing this directly leads to cudaError#11 (Invalid argument)
  struct entry *d_map;
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(&d_map, &d_hm->map, sizeof(d_map), cudaMemcpyDeviceToHost));
  size_t table_size = sizeof(struct entry) * hm->size;
  struct entry *map = (struct entry *) malloc(table_size);
  ASSERT_NO_CUDA_ERROR(
    cudaMemcpy(map, d_map, table_size, cudaMemcpyDeviceToHost));
  hm->map = map;

  return hm;
}

// --------------------------------------------------------------------------
// GPU helper functions
// --------------------------------------------------------------------------

/*
 * Retrieve data from a *device* hashmap within a gpu kernel
 *
 * struct entry *e = cuda_retrieve_hashmap(<hm>, <key>);
 * if (e) { //hit } else { //miss }
 */
__device__ struct entry *cuda_retrieve_hashmap(struct hashmap *hm, int key) {
  struct entry *result = NULL;
  int hash = key % (hm->size);
  if (hm->map[hash].valid && hm->map[hash].key == key) {
    result = &hm->map[hash];
  } else {
    int probe = (hash + 1) % hm->size;
    while (hm->map[probe].key != key && probe != hash) {
      probe = (probe + 1) % hm->size;
    }

    if (probe != hash && hm->map[probe].valid) {
      result = &hm->map[probe];
    }
  }
  return result;
}

// no insertion on device

#endif
