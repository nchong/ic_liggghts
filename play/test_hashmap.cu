#include "hashmap.cu"
#include "cuPrintf.cu"

//print a *device* copy of a hashmap
__global__ void cuda_print_hashmap(struct hashmap *hm) {
  cuPrintf("GPU Hashmap\n");
  for (int i=0; i<hm->size; i++) {
    if (hm->map[i].valid) {
      cuPrintf("%d: %d |-> {%f, %f, %f}\n", 
        i, hm->map[i].key,
        hm->map[i].shear[0],
        hm->map[i].shear[1],
        hm->map[i].shear[2]);
    } else {
      cuPrintf("%d: #\n", i);
    }
  }
}

//dummy update hashmap
__global__ void update_hashmap(struct hashmap *hm) {
  int i = threadIdx.x;
  if (i < hm->size) {
    if (hm->map[i].valid) {
      hm->map[i].shear[0] += 100;
      hm->map[i].shear[1] += 200;
      hm->map[i].shear[2] += 300;
    }
  }
}

int main() {
  struct hashmap *hm = create_hashmap(4);
  double shear0[3] = {1,2,3};
  double shear1[3] = {4,5,6};
  double shear2[3] = {7,8,9};
  double shear3[3] = {10,11,12};
  insert_hashmap(hm, 0, shear0);
  insert_hashmap(hm, 1, shear1);
  insert_hashmap(hm, 2, shear2);
  insert_hashmap(hm, 4, shear3);

  printf("CPU Hashmap(0)\n");
  print_hashmap(hm);

  struct hashmap *d_hm = c_to_device_hashmap(hm);
  memset(hm->map, 0, sizeof(struct entry) * hm->size);
  destroy_hashmap(hm);

  //pre-amble
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("pre-kernel err is %s.\n", cudaGetErrorString(err));
    exit(1);
  }

  cudaPrintfInit();
  update_hashmap<<<1,8>>>(d_hm);
  cuda_print_hashmap<<<1,1>>>(d_hm);
  cudaPrintfDisplay(stdout, true);
  cudaPrintfEnd();

  //post-amble
  cudaThreadSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("post-kernel err is %s.\n", cudaGetErrorString(err));
    exit(1);
  }

  printf("CPU Hashmap(1)\n");
  hm = device_to_c_hashmap(d_hm);
  print_hashmap(hm);

  return 0;
}
