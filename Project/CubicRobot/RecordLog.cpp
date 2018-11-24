#include <iostream>
#include <fstream>
#include "RecordLog.h"
#include "opengl.h"
#include "cube_model.h"
using std::cout;
using std::endl;

std::ofstream ofs_energy("./log/energy.csv");
std::ofstream ofs_time("./log/time.csv");
double total_time=0.0;

void RecordLog(bool verbose){
  if (t==0.0){
      ofs_energy << "time" << "," << "frame" <<","
                 << "Elastic Potential Energy" << ","
                 << "Gravitational Potential Energy" << ","
                 << "Contact Potential Energy" << ","
                 << "Kinetic Energy"  << ","
                 << "External Energy" << ","
                 << "Total Energy" << endl;
  }

  ofs_energy << t << ", " << frame <<", "
             << totalPEk << ", " << totalPEg << ", "
             << totalPEc << ", " << totalKE << ", "
             << totalEE  << ", " << totalEnergy << endl;

  ofs_time << elapsed_seconds.count() << endl;
  total_time += elapsed_seconds.count();
  int num_eval = (int)nt/total_time*N_SPRING;
  if (verbose) cout << "Physics Engine :" << elapsed_seconds.count() << "/ Average" << num_eval << "\n";
}

