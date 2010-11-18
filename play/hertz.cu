#include "cuPrintf.cu"
#include <stdio.h>
#define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875

__global__ void pair_interaction(
  //inputs
    float *xi, float *xj,           //position
    float *vi, float *vj,           //velocity
    float *omegai, float *omegaj,   //rotational velocity
    float radi, float radj,         //radius
    float rmassi, float rmassj,     //] mass
    float massi, float massj,       //]
    int typei, int typej,           //type
    float dt,                       //timestep
   
  //contact model parameters inputs
    int num_atom_types,
    float *Yeff,
    float *Geff,
    float *betaeff,
    float *coeffFrict,
    float nktv2p,

  //inouts
    float *shear,
    float *torque,
    float *force) {

  // del is the vector from j to i
  float delx = xi[0] - xj[0];
  float dely = xi[1] - xj[1];
  float delz = xi[2] - xj[2];

  float rsq = delx*delx + dely*dely + delz*delz;
  float radsum = radi + radj;
  if (rsq >= radsum*radsum) {
    //unset non-touching atoms
  } else {
    //distance between centres of atoms i and j
    //or, magnitude of del vector
    float r = sqrt(rsq); 
    float rinv = 1.0/r;
    float rsqinv = 1.0/rsq;
	
    // relative translational velocity
    float vr1 = vi[0] - vj[0];
    float vr2 = vi[1] - vj[1];
    float vr3 = vi[2] - vj[2];

    // normal component
    float vnnr = vr1*delx + vr2*dely + vr3*delz;
    float vn1 = delx*vnnr * rsqinv;
    float vn2 = dely*vnnr * rsqinv;
    float vn3 = delz*vnnr * rsqinv;

    // tangential component
    float vt1 = vr1 - vn1;
    float vt2 = vr2 - vn2;
    float vt3 = vr3 - vn3;

    // relative rotational velocity
    float wr1 = (radi*omegai[0] + radj*omegaj[0]) * rinv;
    float wr2 = (radi*omegai[1] + radj*omegaj[1]) * rinv;
    float wr3 = (radi*omegai[2] + radj*omegaj[2]) * rinv;

    // normal forces = Hookian contact + normal velocity damping
    float mi,mj;
    if (rmassi && rmassj) {
      mi=rmassi;
      mj=rmassj;
    } else if (massi && massj) {
      mi=massi;
      mj=massj;
    } else {
      //this should never fire
      return;
    }
    float meff = mi*mj/(mi+mj);
    //not-implemented: freeze_group_bit

    float deltan = radsum-r;

    //derive contact model parameters (inlined)
    //Yeff, Geff, betaeff, coeffFrict are lookup tables
    //todo: put these in shared memory
    int typeij = typei + (typej * num_atom_types);
    float reff = radi * radj / (radi + radj);
    float sqrtval = sqrt(reff * deltan);
    float Sn = 2.    * Yeff[typeij] * sqrtval;
    float St = 8.    * Geff[typeij] * sqrtval;
    float kn = 4./3. * Yeff[typeij] * sqrtval;
    float kt = St;
    float gamman=-2.*sqrtFiveOverSix*betaeff[typeij]*sqrt(Sn*meff);
    float gammat=-2.*sqrtFiveOverSix*betaeff[typeij]*sqrt(St*meff);
    float xmu=coeffFrict[typeij];
    kn /= nktv2p;
    kt /= nktv2p;

    //if dampflag gammat = 0
    float damp = gamman*vnnr*rsqinv;  
	  float ccel = kn*(radsum-r)*rinv - damp;

    //not-implemented cohesionflag

    // relative velocities
    float vtr1 = vt1 - (delz*wr2-dely*wr3);
    float vtr2 = vt2 - (delx*wr3-delz*wr1);
    float vtr3 = vt3 - (dely*wr1-delx*wr2);

    // shear history effects
    shear[0] += vtr1 * dt;
    shear[1] += vtr2 * dt;
    shear[2] += vtr3 * dt;

    // rotate shear displacements
    float rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
    rsht *= rsqinv;

    shear[0] -= rsht*delx;
    shear[1] -= rsht*dely;
    shear[2] -= rsht*delz;

    // tangential forces = shear + tangential velocity damping
    float fs1 = - (kt*shear[0] + gammat*vtr1); 
    float fs2 = - (kt*shear[1] + gammat*vtr2); 
    float fs3 = - (kt*shear[2] + gammat*vtr3); 

    // rescale frictional displacements and forces if needed
    float fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    float fn = xmu * fabs(ccel*r);
    if (fs > fn) {
      float shrmag = sqrt(shear[0]*shear[0] + 
                          shear[1]*shear[1] +
                          shear[2]*shear[2]);
      if (shrmag != 0.0) {
        shear[0] = (fn/fs) * (shear[0] + gammat*vtr1/kt) - gammat*vtr1/kt;
        shear[1] = (fn/fs) * (shear[1] + gammat*vtr2/kt) - gammat*vtr2/kt;
        shear[2] = (fn/fs) * (shear[2] + gammat*vtr3/kt) - gammat*vtr3/kt;
        fs1 *= fn/fs;
        fs2 *= fn/fs;
        fs3 *= fn/fs;
      } else {
        fs1 = fs2 = fs3 = 0.0;
      }
    }

    float tor1 = rinv * (dely*fs3 - delz*fs2);
    float tor2 = rinv * (delz*fs1 - delx*fs3);
    float tor3 = rinv * (delx*fs2 - dely*fs1);

    // this is what we've been working up to!
    force[0] += delx*ccel + fs1;
    force[1] += dely*ccel + fs2;
    force[2] += delz*ccel + fs3;

    torque[0] -= radi*tor1;
    torque[1] -= radi*tor2;
    torque[2] -= radi*tor3;
  }
}

int main(void) {
  //-----------
  // Inputs
  //-----------
  float xi[3] = {0.052408, 0.037245, 0.003795};
  float xj[3] = {0.053327,0.035612,0.001654};
  float vi[3] = {0.000000, 0.000000, -1.513257};
  float vj[3] = {0.000000,0.000000,-0.548084};
  float omegai[3] = {0.000000, 0.000000, 0.000000};
  float omegaj[3] = {0.000000,0.000000,0.000000};
  float radi = 0.001036;
  float radj = 0.001815;
  float rmassi = 0.000012;
  float rmassj = 0.000066;
  float massi = 0.0;
  float massj = 0.0;
  int typei = 1;
  int typej = 1;
  float dt = 0.000010;
  //-----------
  int max_type = 2;
  float Yeff[2][2];
  float Geff[2][2];
  float betaeff[2][2];
  float coeffFrict[2][2];
  Yeff[1][1] = 3134796.238245;
  Geff[1][1] = 556173.526140;
  betaeff[1][1] = -0.357857;
  coeffFrict[1][1] = 0.500000;
  float nktv2p = 1.000000;
  //-----------
  // Inouts
  //-----------
  float shear[3] = {0.000000, 0.000000, 0.000000};
  float torque[3] = {0.000000, 0.000000, 0.000000};
  float force[3] = {0.000000, 0.000000, 0.000000};
  //-----------
  // Expected outputs
  //-----------
  float expected_shear[3] = {0.000008, -0.000015, 0.000015};
  float expected_torque[3] = {-0.000014, -0.000008, 0.000000};
  float expected_force[3] = {-0.004426, 0.007860, 0.034589};
  //-----------

  //-----------
  // Device-versions of inputs
  //-----------
  float *d_xi; float *d_xj;
  float *d_vi; float *d_vj;
  float *d_omegai; float *d_omegaj;
  // not required for float inputs (rad, rmass, mass, type and dt)
  //-----------
  float *d_Yeff;
  float *d_Geff;
  float *d_betaeff;
  float *d_coeffFrict;

  cudaMalloc((void**)&d_xi, sizeof(float)*3);
  cudaMalloc((void**)&d_xj, sizeof(float)*3);
  cudaMalloc((void**)&d_vi, sizeof(float)*3);
  cudaMalloc((void**)&d_vj, sizeof(float)*3);
  cudaMalloc((void**)&d_omegai, sizeof(float)*3);
  cudaMalloc((void**)&d_omegaj, sizeof(float)*3);
  cudaMemcpy(d_xi, xi, sizeof(float)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_xj, xj, sizeof(float)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_vi, vi, sizeof(float)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_vj, vj, sizeof(float)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_omegai, omegai, sizeof(float)*3, cudaMemcpyHostToDevice);
  cudaMemcpy(d_omegaj, omegaj, sizeof(float)*3, cudaMemcpyHostToDevice);

  //flatten 2D matrices into contiguous memory
  float *h_Yeff = (float *)malloc(sizeof(float)*max_type*max_type);
  float *h_Geff = (float *)malloc(sizeof(float)*max_type*max_type);
  float *h_betaeff = (float *)malloc(sizeof(float)*max_type*max_type);
  float *h_coeffFrict = (float *)malloc(sizeof(float)*max_type*max_type);
  for (int i=0; i<max_type; i++) {
    for (int j=0; j<max_type; j++) {
      h_Yeff[i + (j*max_type)] = Yeff[i][j];
      h_Geff[i + (j*max_type)] = Geff[i][j];
      h_betaeff[i + (j*max_type)] = betaeff[i][j];
      h_coeffFrict[i + (j*max_type)] = coeffFrict[i][j];
    }
  }
  cudaMalloc((void**)&d_Yeff, sizeof(float)*max_type*max_type);
  cudaMalloc((void**)&d_Geff, sizeof(float)*max_type*max_type);
  cudaMalloc((void**)&d_betaeff, sizeof(float)*max_type*max_type);
  cudaMalloc((void**)&d_coeffFrict, sizeof(float)*max_type*max_type);
  cudaMemcpy(d_Yeff, h_Yeff, sizeof(float)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Geff, h_Geff, sizeof(float)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_betaeff, h_betaeff, sizeof(float)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_coeffFrict, h_coeffFrict, sizeof(float)*max_type*max_type, cudaMemcpyHostToDevice);

  //inouts
  float *d_shear;
  cudaMalloc((void**)&d_shear, sizeof(float)*3);
	cudaMemcpy(d_shear, shear, sizeof(float)*3, cudaMemcpyHostToDevice);

  float *d_torque;
  cudaMalloc((void**)&d_torque, sizeof(float)*3);
	cudaMemcpy(d_torque, torque, sizeof(float)*3, cudaMemcpyHostToDevice);

  float *d_force;
  cudaMalloc((void**)&d_force, sizeof(float)*3);
	cudaMemcpy(d_force, force, sizeof(float)*3, cudaMemcpyHostToDevice);

  cudaError_t err = cudaGetLastError();
  printf("pre-kernel err is %s.\n", cudaGetErrorString(err));

  cudaPrintfInit();

  pair_interaction<<<1,1>>>(
    d_xi,d_xj,
    d_vi,d_vj,
    d_omegai, d_omegaj,
    radi, radj,
    rmassi, rmassj,
    massi, massj,
    typei, typej,
    dt,

    max_type,
    d_Yeff,
    d_Geff,
    d_betaeff,
    d_coeffFrict,
    nktv2p,

    d_shear, d_torque, d_force);

  cudaPrintfDisplay(stdout, true);
  cudaPrintfEnd();

  cudaThreadSynchronize();
  err = cudaGetLastError();
  printf("kernel run was %s.\n", cudaGetErrorString(err));

  cudaMemcpy(shear, d_shear, sizeof(float)*3, cudaMemcpyDeviceToHost);
  cudaMemcpy(torque, d_torque, sizeof(float)*3, cudaMemcpyDeviceToHost);
  cudaMemcpy(force, d_force, sizeof(float)*3, cudaMemcpyDeviceToHost);

  printf("expected shear    = {%f, %f, %f}\n",
    expected_shear[0], expected_shear[1], expected_shear[2]);
  printf("calculated shear  = {%f, %f, %f}\n",
    shear[0], shear[1], shear[2]);

  printf("expected torque   = {%f, %f, %f}\n",
    expected_torque[0], expected_torque[1], expected_torque[2]);
  printf("calculated torque = {%f, %f, %f}\n",
    torque[0], torque[1], torque[2]);

  printf("expected force    = {%f, %f, %f}\n",
    expected_force[0], expected_force[1], expected_force[2]);
  printf("calculated force  = {%f, %f, %f}\n",
    force[0], force[1], force[2]);

  return 0;
}
