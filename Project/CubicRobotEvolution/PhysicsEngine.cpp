#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "CubeGenerator.h"
#include "OpenGL.h"
#include "Vector3D.h"
#include "utility.h"

using std::vector;
using std::cout;
using std::endl;

int nt = 0;
double t = 0.0;
double dt = 0.0005;

// Energy
double totalPEk;
double totalPEc;
double totalPEg;
double totalKE;
double totalEE;
double totalEnergy;
double elapsed_seconds;

// Gravitational acceleration
const double ga = 9.81;

// Contact parameter
const double contactK = 1e5;
const double FrictionCoefficient = 1.0;
vector<bool> slip(N_MASS, true);

// Damping
const double totalDamping = 0.999;

/*
const double springDampingC = 1e1;
const double dragConstant = -1e-2
const double contactDampingC = 1e-1;
*/

// -----------------------------------------------------------
// Physics Engine
// -----------------------------------------------------------

void PhysicsEngine(void){
  Timer tm;

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

  /*
  // Switch damping ON/OFF
  double springC; double contactC; double allC;
  if (!_Damping){
    springC=0.0; contactC=0.0; allC=0.0;
  }else{
    springC=springDampingC; contactC=contactDampingC;
  }
  */
  double allC;
  if (!_Damping){
    allC=0.0;
  }else{
    allC=totalDamping;
  }


  // Spring Force ---------------------------------------------------------------
  for(int i=0; i<N_SPRING; i++){

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
    force[index1] = add(force[index1], spring_force);
    force[index2] = sub(force[index2], spring_force);

    // Damping term
    /*
    vector<double> rel_velo = sub(mass[index2].v, mass[index1].v);
    double damping_force_abs = springC * dotProduct(rel_velo, spring_unit_vec);
    vector<double> damping_force = scaling(spring_unit_vec, damping_force_abs);
    force[index1] = add(force[index1], damping_force);
    force[index2] = sub(force[index2], damping_force);
    */

    totalPEk += 0.5 * spring[i].k * pow((length - spring[i].l0), 2);
  }

  //#pragma omp parallel for
  for(int i=0; i<N_MASS; i++){

    // External Force --------------------------------------------------------------
    if (i==nudgeID){
      force[i] = add(force[i], nudgeForce);
      nudgeForce = {0.0, 0.0, 0.0}; // Reset force.
    }
    //totalEE += 0.0;

    // Gravitational Force -------------------------------------------------------
    vector<double>grav_force = {0.0, 0.0, -mass[i].m*ga};
    force[i] = add(force[i], grav_force);
    totalPEg += mass[i].m * ga * mass[i].p[2];

    // Drag Force ----------------------------------------------------------------
    /*
    if (_Damping){
      vector<double> drag_force = scaling(mass[i].v, dragConstant);
      force[i] = add(force[i], drag_force);
    }
    */

    // Contact Force --------------------------------------------------------------
    if(mass[i].p[2] < 0.0){ //nodeRadius

      // Normal Force
      double normal_force_abs =
            contactK * pow((nodeRadius - mass[i].p[2]), 2);
      /*
      double normal_force_abs =
            contactK * pow((nodeRadius - mass[i].p[2]), 2) -
            contactC * mass[i].v[2];
      */

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
      totalPEc += 1.0/3.0 * contactK * pow((- mass[i].p[2]), 3);

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
    totalKE += 0.5 * mass[i].m * calcNorm(mass[i].v) * calcNorm(mass[i].v);
  }

  totalEnergy = totalPEk + totalPEg + totalPEc + totalKE + totalEE;

  elapsed_seconds =   tm.elapsed();

}


void EnergyRecord(void){
  static std::ofstream ofs_energy("./log/energy.csv");
  static std::ofstream ofs_time("./log/time.csv");

  if (nt==0){
    ofs_energy << "time" << "," << "frame" <<","
               << "Elastic Potential Energy" << ","
               << "Gravitational Potential Energy" << ","
               << "Contact Potential Energy" << ","
               << "Kinetic Energy"  << ","
               << "External Energy" << ","
               << "Total Energy" << endl;

    ofs_time << "nt, Elapsed time, Num. of Spring Evaluations" << endl;
  }


  ofs_energy << t << ", " << frame <<", "
             << totalPEk << ", " << totalPEg << ", "
             << totalPEc << ", " << totalKE << ", "
             << totalEE  << ", " << totalEnergy << endl;

  ofs_time << std::setw(6) << std::right << nt << ", ";
  ofs_time << std::setprecision(6) << elapsed_seconds << ", ";
  ofs_time << std::setw(8) << std::right << (int) (N_SPRING / elapsed_seconds) << endl;

  //total_time += elapsed_seconds;
  //double average_time = total_time / (double)(nt+1);
  //int num_eval = (int) (N_SPRING / average_time);
  //if (verbose) cout << "Physics Engine :" << elapsed_seconds.count() << "/ Average" << num_eval << "\n";

}

void TrajectoryRecord(int slot){

  string str;

  //if(nt==0){
  //  for(int i=0; i<10; i++){
  //    str = "./log/robots/robot" + std::to_string(i) + ".csv";
  //    std::ifstream fs(str);
  //    if (!fs.is_open()) break;
  //  }
  //}

  if (nt==0) str = "./log/robots/robot" + std::to_string(slot) + ".csv";
  static std::ofstream ofs_trajectory(str);

  if(nt==0){
    //cout << "A file: "<< str << " is created." << endl;
    ofs_trajectory << "time" << ","
                   << "id" << ","
                   << "x" << ","
                   << "y" << ","
                   << "z" << endl;
  }

  for (int i=0; i<N_MASS; i++){
    ofs_trajectory << t << ","
                   << i << ","
                   << mass[i].p[0] << ","
                   << mass[i].p[1] << ","
                   << mass[i].p[2] << endl;
  }

}
