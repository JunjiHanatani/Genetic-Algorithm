#include <math.h>
#include <iostream>
#include <omp.h>
#include <chrono>
#include "CubeGenerator.h"
#include "OpenGL.h"
#include "Vector3D.h"
#include "RecordLog.h"
#include "utility.h"

using std::vector;
using std::cout;
using std::endl;

std::chrono::duration<double> elapsed_seconds;

// Energy
double totalPEk;
double totalPEc;
double totalPEg;
double totalKE;
double totalEE;
double totalEnergy;

// Gravitational acceleration
const double ga = 9.81;

// Contact parameter
const double contactK = 1e5;
const double FrictionCoefficient = 0.5;
vector<bool> slip(N_MASS, true);

// Damping
const double totalDamping = 0.999;

// -----------------------------------------------------------
// Physics Engine
// -----------------------------------------------------------

void PhysicsEngine(void){
  auto time0 = std::chrono::system_clock::now();

  // Initialization ------------------------------------------------------------
  vector<vector<double>> force(N_MASS, vector<double>(3, 0.0));
  slip.assign(N_MASS, true);

  // Initialize Energy
  totalPEk=0.0;
  totalPEg=0.0;
  totalPEc=0.0;
  totalKE=0.0;
  double localPEk=0.0;
  double localPEg=0.0;
  double localPEc=0.0;
  double localKE=0.0;

  // Switch friction ON/OFF
  double mu;
  if (!_Friction){
    mu = 0.0;
  }else{
    mu = FrictionCoefficient;
  }

  // Switch damping ON/OFF
  double allC;
  if (!_Damping){
    allC=0.0;
  }else{
    allC=totalDamping;
  }

  // Spring Force ---------------------------------------------------------------
  double local_vec[40*3];
  double global_vec[40*3*4];

  for(int iter=0; iter<40; iter++){

    int i = rank*40+iter;

    if (i<N_SPRING){

      if(_Breathe){
        if (t<0){
          spring[i].l0 = springInitialRestlength[i] +
                         ( breathe_offset[i] + breathe_amp[i] * sin(breathe_phase[i]) ) * (1.0 + t);
        }else{
          spring[i].l0 = springInitialRestlength[i] + breathe_offset[i] +
                         breathe_amp[i] * sin(breathe_omega*t + breathe_phase[i]);
        }
      }

      // A pair of masses associated with the i-th spring.
      int index1 = spring[i].masses[0];
      int index2 = spring[i].masses[1];

      // Spring term
      double length = calcDistance(mass[index1].p, mass[index2].p);
      double spring_force_abs = spring[i].k * (length - spring[i].l0);
      vector<double> spring_unit_vec = calcUnitVector2(mass[index1].p, mass[index2].p);
      vector<double> spring_force = scaling(spring_unit_vec, spring_force_abs);

      local_vec[iter] = spring_force[0];
      local_vec[iter+1] = spring_force[1];
      local_vec[iter+2] = spring_force[2];

      localPEk += 0.5 * spring[i].k * pow((length - spring[i].l0), 2);
    }

  }

  MPI_Gather( &local_vec, 40*3, MPI_DOUBLE, &global_vec, 40*3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Reduce( &localPEk, &totalPEk, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

  for(int i=0; i<N_SPRING; i++){
    int index1 = spring[i].masses[0];
    int index2 = spring[i].masses[1];
    for (int j=0; j<3; j++){
      force[index1][j] += global_vec[i*3 + j];
      force[index2][j] -= global_vec[i*3 + j];
    }
  }

/* MASS */

  for(int iter=0; iter<16; iter++){

    int i = rank*16+iter;

    if (i<N_MASS){

      // Gravitational Force -------------------------------------------------------
      vector<double>grav_force = {0.0, 0.0, -mass[i].m*ga};
      force[i] = add(force[i], grav_force);
      localPEg += mass[i].m * ga * mass[i].p[2];

      // Contact Force --------------------------------------------------------------
      if(mass[i].p[2] < 0.0){ //nodeRadius

        // Normal Force
        double normal_force_abs =
              contactK * pow((nodeRadius - mass[i].p[2]), 2);

        vector<double> normal_force = {0.0, 0.0, normal_force_abs};

        // Friction Force
        vector<double> horizontal_force = {force[i][0], force[i][1], 0.0};
        double horizontal_force_abs = calcNorm(horizontal_force);

        vector<double> friction_force;
        double friction_force_abs = normal_force_abs * mu;

        if (horizontal_force_abs > friction_force_abs){ // Slip
          slip[i] = true;
          friction_force = calcVector(horizontal_force, -friction_force_abs);
        }else{ // Stick
          slip[i] = false;
          friction_force = scaling(horizontal_force, -1.0);
        }

        // Contact Force (Normal Force + Friction Force)
        vector<double> contact_force = add(normal_force, friction_force);
        force[i] = add(force[i], contact_force);

        // Contact Potential Energy
        localPEc += 1.0/3.0 * contactK * pow((nodeRadius - mass[i].p[2]), 3);

      }

    // Solve equations of motion.
    for (int j=0; j<3; j++){
      mass[i].a[j] = force[i][j]/mass[i].m;
      mass[i].v[j] += mass[i].a[j] * dt;
      if (_Damping) mass[i].v[j] = allC * mass[i].v[j];
      if (!slip[i] && (j==0 || j==1)) mass[i].v[j]=0.0;
      mass[i].p[j] += mass[i].v[j] * dt;
    }

    // Kinetic Energy
    localKE += 0.5 * mass[i].m * calcNorm(mass[i].v) * calcNorm(mass[i].v);
  }

  totalEnergy = totalPEk + totalPEg + totalPEc + totalKE;

  auto time1 = std::chrono::system_clock::now();
  elapsed_seconds = time1-time0;

}

