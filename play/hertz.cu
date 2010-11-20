#include "cuPrintf.cu"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

// --------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------

// testcase datastructure
struct params {
  double xi[3]; double xj[3];
  double vi[3]; double vj[3];
  double omegai[3]; double omegaj[3];
  double radi; double radj;
  double rmassi; double rmassj;
  double massi; double massj;
  double shear[3]; double torque[3]; double force[3];
  double expected_shear[3]; double expected_torque[3]; double expected_force[3];
};

struct params parse_csv_string(const char *str) {
  struct params result;
  sscanf(str,
    "%lf, %lf, %lf, %lf, %lf, %lf, "                   //x
    "%lf, %lf, %lf, %lf, %lf, %lf, "                   //v
    "%lf, %lf, %lf, %lf, %lf, %lf, "                   //omega
    "%lf, %lf, "                                       //radius
    "%lf, %lf, "                                       //rmass
    "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, "    //shear, torque, force
    "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, ",   //*expected* shear, torque, force
    &result.xi[0], &result.xi[1], &result.xi[2],
    &result.xj[0], &result.xj[1], &result.xj[2],
    &result.vi[0], &result.vi[1], &result.vi[2],
    &result.vj[0], &result.vj[1], &result.vj[2],
    &result.omegai[0], &result.omegai[1], &result.omegai[2],
    &result.omegaj[0], &result.omegaj[1], &result.omegaj[2],
    &result.radi, &result.radj,
    &result.rmassi, &result.rmassj,
    &result.shear[0], &result.shear[1], &result.shear[2],
    &result.torque[0], &result.torque[1], &result.torque[2],
    &result.force[0], &result.force[1], &result.force[2],
    &result.expected_shear[0], &result.expected_shear[1], &result.expected_shear[2],
    &result.expected_torque[0], &result.expected_torque[1], &result.expected_torque[2],
    &result.expected_force[0], &result.expected_force[1], &result.expected_force[2]
    );
  result.massi = result.massj = 0;
  return result;
}

bool check_result_vector(const char* id, double expected[3], double actual[3], const double epsilon) {
  static bool verbose = true;
  bool flag = (fabs(expected[0] - actual[0]) > epsilon ||
               fabs(expected[1] - actual[1]) > epsilon ||
               fabs(expected[2] - actual[2]) > epsilon);
  const char *marker = flag ? "***" : "   ";

  if (flag || verbose) {
    printf("%s%s: {%.16f, %.16f, %.16f} / {%.16f, %.16f, %.16f}%s\n",
        marker,
        id,
        expected[0], expected[1], expected[2],
        actual[0], actual[1], actual[2],
        marker
        );
  }
  return flag;
}

// --------------------------------------------------------------------------
// GPU Kernel
// --------------------------------------------------------------------------
#define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875

__global__ void pair_interaction(
  //inputs
    double *xi, double *xj,           //position
    double *vi, double *vj,           //velocity
    double *omegai, double *omegaj,   //rotational velocity
    double radi, double radj,         //radius
    double rmassi, double rmassj,     //] mass
    double massi, double massj,       //]
    int typei, int typej,           //type
    double dt,                       //timestep

  //contact model parameters inputs
    int num_atom_types,
    double *Yeff,
    double *Geff,
    double *betaeff,
    double *coeffFrict,
    double nktv2p,

  //inouts
    double *shear,
    double *torque,
    double *force) {

  // del is the vector from j to i
  double delx = xi[0] - xj[0];
  double dely = xi[1] - xj[1];
  double delz = xi[2] - xj[2];

  double rsq = delx*delx + dely*dely + delz*delz;
  double radsum = radi + radj;
  if (rsq >= radsum*radsum) {
    //unset non-touching atoms
    shear[0] = 0.0;
    shear[1] = 0.0;
    shear[2] = 0.0;
  } else {
    //distance between centres of atoms i and j
    //or, magnitude of del vector
    double r = sqrt(rsq);
    double rinv = 1.0/r;
    double rsqinv = 1.0/rsq;
	
    // relative translational velocity
    double vr1 = vi[0] - vj[0];
    double vr2 = vi[1] - vj[1];
    double vr3 = vi[2] - vj[2];

    // normal component
    double vnnr = vr1*delx + vr2*dely + vr3*delz;
    double vn1 = delx*vnnr * rsqinv;
    double vn2 = dely*vnnr * rsqinv;
    double vn3 = delz*vnnr * rsqinv;

    // tangential component
    double vt1 = vr1 - vn1;
    double vt2 = vr2 - vn2;
    double vt3 = vr3 - vn3;

    // relative rotational velocity
    double wr1 = (radi*omegai[0] + radj*omegaj[0]) * rinv;
    double wr2 = (radi*omegai[1] + radj*omegaj[1]) * rinv;
    double wr3 = (radi*omegai[2] + radj*omegaj[2]) * rinv;

    // normal forces = Hookian contact + normal velocity damping
    double mi,mj;
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
    double meff = mi*mj/(mi+mj);
    //not-implemented: freeze_group_bit

    double deltan = radsum-r;

    //derive contact model parameters (inlined)
    //Yeff, Geff, betaeff, coeffFrict are lookup tables
    //todo: put these in shared memory
    int typeij = typei + (typej * num_atom_types);
    double reff = radi * radj / (radi + radj);
    double sqrtval = sqrt(reff * deltan);
    double Sn = 2.    * Yeff[typeij] * sqrtval;
    double St = 8.    * Geff[typeij] * sqrtval;
    double kn = 4./3. * Yeff[typeij] * sqrtval;
    double kt = St;
    double gamman=-2.*sqrtFiveOverSix*betaeff[typeij]*sqrt(Sn*meff);
    double gammat=-2.*sqrtFiveOverSix*betaeff[typeij]*sqrt(St*meff);
    double xmu=coeffFrict[typeij];
    //not-implemented if (dampflag == 0) gammat = 0;
    kn /= nktv2p;
    kt /= nktv2p;

    double damp = gamman*vnnr*rsqinv;
	  double ccel = kn*(radsum-r)*rinv - damp;

    //not-implemented cohesionflag

    // relative velocities
    double vtr1 = vt1 - (delz*wr2-dely*wr3);
    double vtr2 = vt2 - (delx*wr3-delz*wr1);
    double vtr3 = vt3 - (dely*wr1-delx*wr2);

    // shear history effects
    shear[0] += vtr1 * dt;
    shear[1] += vtr2 * dt;
    shear[2] += vtr3 * dt;

    // rotate shear displacements
    double rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
    rsht *= rsqinv;

    shear[0] -= rsht*delx;
    shear[1] -= rsht*dely;
    shear[2] -= rsht*delz;

    // tangential forces = shear + tangential velocity damping
    double fs1 = - (kt*shear[0] + gammat*vtr1);
    double fs2 = - (kt*shear[1] + gammat*vtr2);
    double fs3 = - (kt*shear[2] + gammat*vtr3);

    // rescale frictional displacements and forces if needed
    double fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    double fn = xmu * fabs(ccel*r);
    double shrmag = 0;
    if (fs > fn) {
      shrmag = sqrt(shear[0]*shear[0] +
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

    double fx = delx*ccel + fs1;
    double fy = dely*ccel + fs2;
    double fz = delz*ccel + fs3;

    double tor1 = rinv * (dely*fs3 - delz*fs2);
    double tor2 = rinv * (delz*fs1 - delx*fs3);
    double tor3 = rinv * (delx*fs2 - dely*fs1);

    // this is what we've been working up to!
    force[0] += fx;
    force[1] += fy;
    force[2] += fz;

    torque[0] -= radi*tor1;
    torque[1] -= radi*tor2;
    torque[2] -= radi*tor3;
  }
}

// --------------------------------------------------------------------------
// Main
// --------------------------------------------------------------------------
int main(void) {
  // Stiffness lookup tables (indexed on atom type)
  double Yeff[2][2];
  double Geff[2][2];
  double betaeff[2][2];
  double coeffFrict[2][2];

  // Inputs fixed across testcases
  double dt = 0.000010;
  Yeff[1][1] = 3134796.2382445144467056;
  Geff[1][1] = 556173.5261401557363570;
  betaeff[1][1] = -0.3578571305033167;
  coeffFrict[1][1] = 0.5000000000000000;
  double nktv2p = 1.000000;
  int typei = 1;
  int typej = 1;
  int max_type = 2;

  // Device versions of inputs
  double *d_xi; double *d_xj;
  double *d_vi; double *d_vj;
  double *d_omegai; double *d_omegaj;
  // not required for double inputs (rad, rmass, mass, type and dt)
  double *d_Yeff;
  double *d_Geff;
  double *d_betaeff;
  double *d_coeffFrict;
  double *d_shear;
  double *d_torque;
  double *d_force;

  // Device allocation of input parameters
  cudaMalloc((void**)&d_xi, sizeof(double)*3);
  cudaMalloc((void**)&d_xj, sizeof(double)*3);
  cudaMalloc((void**)&d_vi, sizeof(double)*3);
  cudaMalloc((void**)&d_vj, sizeof(double)*3);
  cudaMalloc((void**)&d_omegai, sizeof(double)*3);
  cudaMalloc((void**)&d_omegaj, sizeof(double)*3);
  cudaMalloc((void**)&d_Yeff, sizeof(double)*max_type*max_type);
  cudaMalloc((void**)&d_Geff, sizeof(double)*max_type*max_type);
  cudaMalloc((void**)&d_betaeff, sizeof(double)*max_type*max_type);
  cudaMalloc((void**)&d_coeffFrict, sizeof(double)*max_type*max_type);
  cudaMalloc((void**)&d_shear, sizeof(double)*3);
  cudaMalloc((void**)&d_torque, sizeof(double)*3);
  cudaMalloc((void**)&d_force, sizeof(double)*3);

  // Flatten 2D lookup tables into contiguous memory
  double *h_Yeff = (double *)malloc(sizeof(double)*max_type*max_type);
  double *h_Geff = (double *)malloc(sizeof(double)*max_type*max_type);
  double *h_betaeff = (double *)malloc(sizeof(double)*max_type*max_type);
  double *h_coeffFrict = (double *)malloc(sizeof(double)*max_type*max_type);
  for (int i=0; i<max_type; i++) {
    for (int j=0; j<max_type; j++) {
      h_Yeff[i + (j*max_type)] = Yeff[i][j];
      h_Geff[i + (j*max_type)] = Geff[i][j];
      h_betaeff[i + (j*max_type)] = betaeff[i][j];
      h_coeffFrict[i + (j*max_type)] = coeffFrict[i][j];
    }
  }
  cudaMemcpy(d_Yeff, h_Yeff, sizeof(double)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Geff, h_Geff, sizeof(double)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_betaeff, h_betaeff, sizeof(double)*max_type*max_type, cudaMemcpyHostToDevice);
  cudaMemcpy(d_coeffFrict, h_coeffFrict, sizeof(double)*max_type*max_type, cudaMemcpyHostToDevice);

  // Open testcase datafile
  std::ifstream data("pairwise_data.csv", std::fstream::in);
  if (!data.is_open()) {
    printf("Could not find/open file [pairwise_data.csv]\n");
    exit(-1);
  }

  // Test loop over datafile
  std::string line;
  while(std::getline(data, line)) {
    struct params testcase = parse_csv_string(line.c_str());
    cudaMemcpy(d_xi, testcase.xi, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_xj, testcase.xj, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vi, testcase.vi, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vj, testcase.vj, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_omegai, testcase.omegai, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_omegaj, testcase.omegaj, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_shear, testcase.shear, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_torque, testcase.torque, sizeof(double)*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_force, testcase.force, sizeof(double)*3, cudaMemcpyHostToDevice);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      printf("pre-kernel err is %s.\n", cudaGetErrorString(err));
      exit(-1);
    }

    cudaPrintfInit();
    pair_interaction<<<1,1>>>(
      d_xi,d_xj,
      d_vi,d_vj,
      d_omegai, d_omegaj,
      testcase.radi, testcase.radj,
      testcase.rmassi, testcase.rmassj,
      testcase.massi, testcase.massj,
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
    if (err != cudaSuccess) {
      printf("post-kernel err is %s.\n", cudaGetErrorString(err));
      exit(-1);
    }

    // Check results
    double shear[3];
    double torque[3];
    double force[3];
    cudaMemcpy(shear, d_shear, sizeof(double)*3, cudaMemcpyDeviceToHost);
    cudaMemcpy(torque, d_torque, sizeof(double)*3, cudaMemcpyDeviceToHost);
    cudaMemcpy(force, d_force, sizeof(double)*3, cudaMemcpyDeviceToHost);

    const double epsilon = 0.00001;
    check_result_vector("shear ", testcase.expected_shear, shear, epsilon);
    check_result_vector("torque", testcase.expected_torque, torque, epsilon);
    check_result_vector("force ", testcase.expected_force, force, epsilon);
  }

  return 0;
}
