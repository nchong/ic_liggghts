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
  hm->size = size;
  hm->map = (struct entry *) malloc(sizeof(struct entry) * size);
  for (int i=0; i<size; i++) {
    hm->map[i].valid = false;
    hm->map[i].key = 0;
    hm->map[i].shear[0] = 0;
    hm->map[i].shear[1] = 0;
    hm->map[i].shear[2] = 0;
  }
  return hm;
}

void destroy_hashmap(struct hashmap *hm) {
  free(hm->map);
  free(hm);
}

struct entry *insert_hashmap(struct hashmap *hm, int key, double shear[3]) {
  int hash = key % (hm->size);

  //linear probe on collision
  if (hm->map[hash].valid) {
    int probe = hash + 1;
    while (hm->map[probe].valid && probe != hash) {
      probe = (probe + 1) % hm->size;
    }
    
    if (probe == hash) {
      printf("hm is full!\n"); //could resize here
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
  struct entry *result = NULL;
  int hash = key % (hm->size);
  if (hm->map[hash].valid && hm->map[hash].key == key) {
    result = &hm->map[hash];
  } else {
    int probe = hash + 1;
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
  cudaMalloc((void **)&d_map, table_size);
  cudaMemcpy(d_map, hm->map, table_size, cudaMemcpyHostToDevice);

  struct hashmap *d_hm;
  cudaMalloc((void **)&d_hm, sizeof(struct hashmap));
  cudaMemcpy(&d_hm->size, &hm->size, sizeof(int), cudaMemcpyHostToDevice);
  //nb copying sizeof pointer *not* the table here
  cudaMemcpy(&d_hm->map, &d_map, sizeof(&d_map), cudaMemcpyHostToDevice);

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("create_device_hashmap_of: %s.\n", cudaGetErrorString(err));
    exit(1);
  }

  return d_hm;
}

/*
 * Inverse of c_to_device_hashmap. Create a C copy of a device hashmap
 */
struct hashmap *device_to_c_hashmap(struct hashmap *d_hm) {
  struct hashmap *hm = (struct hashmap *) malloc(sizeof(struct hashmap));
  cudaMemcpy(&hm->size, &d_hm->size, sizeof(int), cudaMemcpyDeviceToHost);

  //nb this is a *device* pointer to unwrap the hashmap struct
  //doing this directly leads to cudaError#11 (Invalid argument)
  struct entry *d_map;
  cudaMemcpy(&d_map, &d_hm->map, sizeof(d_map), cudaMemcpyDeviceToHost);
  size_t table_size = sizeof(struct entry) * hm->size;
  struct entry *map = (struct entry *) malloc(table_size);
  cudaMemcpy(map, d_map, table_size, cudaMemcpyDeviceToHost);
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
    int probe = hash + 1;
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
