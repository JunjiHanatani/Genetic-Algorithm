#include <math.h>
#include <iostream>
#include <omp.h>
#include <chrono>
#include "cube_model.h"
#include "opengl.h"
#include "Vector3D.h"
#include "RecordLog.h"

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

// Structure
const int N_MASS=8;
const int N_SPRING=28;
Mass mass[N_MASS];
Spring spring[N_SPRING];

// Node parameter
const double nodeRadius = 0.01;
const double nodeMass = 0.1;
vector<double> nodeInitialPosition = {0.0, 0.0, nodeRadius+0.01};

// Edge parameter
const double springK = 1e3;
const double springDampingC = 1e-1;
vector<double> initial_length(N_SPRING);

// Contact parameter
const double contactK = 1e5;
const double contactDampingC = 1e1;
const double FrictionCoefficient = 0.5;
vector<bool> slip(N_MASS, true);

// Actuator
const double breathe_amp = 0.03;
const double breathe_freq = 1.0;

void InitializeCube(void){

  mass[0].p = {0.0, 0.0, 0.0};
  mass[1].p = {0.1, 0.0, 0.0};
  mass[2].p = {0.1, 0.1, 0.0};
  mass[3].p = {0.0, 0.1, 0.0};
  mass[4].p = {0.0, 0.0, 0.1};
  mass[5].p = {0.1, 0.0, 0.1};
  mass[6].p = {0.1, 0.1, 0.1};
  mass[7].p = {0.0, 0.1, 0.1};

    // Set initial mass conditions.
  for (int i=0; i<N_MASS; i++){
      mass[i].m = nodeMass;
      mass[i].p = add(mass[i].p, nodeInitialPosition);
      mass[i].v = {0.0, 0.0, 0.0};
      mass[i].a = {0.0, 0.0, 0.0};
  }

  // Set initial spring conditions.
  int k=0;
  for (int i=0; i<N_MASS; i++){
    for (int j=i+1; j<N_MASS; j++){
      spring[k].masses = {i, j};
      spring[k].l0 = calcDistance(mass[i].p, mass[j].p);
      spring[k].k = springK;
      k += 1;
    }
  }

  //  Store initial spring lengths.
  for (int i=0; i<N_SPRING; i++){
    initial_length[i] = spring[i].l0;
  }

}

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
  totalEE=0.0;

  // Switch friction ON/OFF
  double mu;
  if (!_Friction){
    mu = 0.0;
  }else{
    mu = FrictionCoefficient;
  }

  // Switch damping ON/OFF
  double springC; double contactC;
  if (!_Damping){
    springC=0.0; contactC=0.0;
  }else{
    springC=springDampingC; contactC=contactDampingC;
  }

  // Spring Force ---------------------------------------------------------------
  for(int i=0; i<N_SPRING; i++){

    if (_Breathe){
        spring[i].l0 = initial_length[i] + breathe_amp * sin(2.0*PI*breathe_freq*t);
    }

    // A pair of masses associated with the i-th spring.
    int index1 = spring[i].masses[0];
    int index2 = spring[i].masses[1];

    // Spring term
    double length = calcDistance(mass[index1].p, mass[index2].p);
    double spring_force_abs = spring[i].k * (length - spring[i].l0);
    vector<double> spring_unit_vec = calcUnitVector2(mass[index1].p, mass[index2].p);
    vector<double> spring_force = scaling(spring_unit_vec, spring_force_abs);
    force[index1] = add(force[index1], spring_force);
    force[index2] = sub(force[index2], spring_force);

    // Damping term
    vector<double> rel_velo = sub(mass[index2].v, mass[index1].v);
    double damping_force_abs = springC * dotProduct(rel_velo, spring_unit_vec);
    vector<double> damping_force = scaling(spring_unit_vec, damping_force_abs);
    force[index1] = add(force[index1], damping_force);
    force[index2] = sub(force[index2], damping_force);
    totalPEk += 0.5 * spring[i].k * pow((length - spring[i].l0), 2);
  }


  //#pragma omp parallel for
  for(int i=0; i<N_MASS; i++){

    // External Force --------------------------------------------------------------
    if (i==nudgeID){
      force[i] = add(force[i], nudgeForce);
      nudgeForce = {0.0, 0.0, 0.0}; // Reset force.
    }
    totalEE += 0.0;

    // Gravitational Force -------------------------------------------------------
    vector<double>grav_force = {0.0, 0.0, -mass[i].m*ga};
    force[i] = add(force[i], grav_force);
    totalPEg += mass[i].m * ga * mass[i].p[2];

    // Contact Force --------------------------------------------------------------
    if(mass[i].p[2] < nodeRadius){

      // Normal Force
      double normal_force_abs =
            contactK * pow((nodeRadius - mass[i].p[2]), 2) -
            contactC * mass[i].v[2];

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
      totalPEc += 1.0/3.0 * contactK * pow((nodeRadius - mass[i].p[2]), 3);

    }

    // Solve equations of motion.
    for (int j=0; j<3; j++){
      mass[i].a[j] = force[i][j]/mass[i].m;
      mass[i].v[j] += mass[i].a[j] * dt;
      if (!slip[i]) {mass[i].v[0]=0.0; mass[i].v[1]=0.0;}
      mass[i].p[j] += mass[i].v[j] * dt;
    }

    // Kinetic Energy
    totalKE += 0.5 * mass[i].m * calcNorm(mass[i].v) * calcNorm(mass[i].v);
  }

  totalEnergy = totalPEk + totalPEg + totalPEc + totalKE + totalEE;

  auto time1 = std::chrono::system_clock::now();
  elapsed_seconds = time1-time0;

  RecordLog(false);

}
