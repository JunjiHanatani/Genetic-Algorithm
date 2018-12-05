#include <iostream>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <vector>
#include <math.h>
#include <string>
#include "Vector3D.h"
#include "CubeGenerator.h"
#include "utility.h"
#include "test.h"
#include "GeneticOperators.h"
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::vector;

void test(void){

  std::ofstream ofs_test("./test.csv");

  // Define the covariance matrix and the mean

  const int N = 5;
  Eigen::Matrix3d sigma[N];
  Eigen::Vector3d mean[N];
  double thresh = 0.0;

  for(int i=0; i<N; i++){
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(3,3) * sqrt(robot_diameter/6.0);
    sigma[i] = m.transpose() * m;
    mean[i] << get_rand_range_dbl(range_min[0], range_max[0]),
               get_rand_range_dbl(range_min[1], range_max[1]),
               get_rand_range_dbl(range_min[2], range_max[2]);
  }


  Mvn mvn[N];
  for (int i=0; i<N; i++){
    mvn[i].Set(mean[i], sigma[i]);
  }

  for(Mass m: mass){
    Eigen::Vector3d pt;
    pt << m.p[0], m.p[1], m.p[2];
    double w = 0.0;
    for (int l=0; l<N; l++) w += mvn[l].pdf(pt);
    w = w/(double)N;
    if(w > thresh) ofs_test << pt(0) <<"," << pt(1) << "," << pt(2) << "," << w << endl;
  }

}
