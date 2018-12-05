#include <iostream>
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


using std::cout;
using std::endl;
using std::vector;

void test(void){

  double x = 1.0;
  int num[5] = {0, 1, 2, 3, 4};
  double* p = &x;
  int* pi;

  cout << num << endl;
  pi = num;

  cout << p << endl;
  cout << pi[0] << endl;

  for(int& i:num) i = i + 1;
  //for (int i=0; i<5; i++){
    //num[i] += 1;
    //int& j = num[i];
    //j += 1;
  //}

  cout << num[0] << endl;
  cout << num[1] << endl;
  cout << num[2] << endl;

}

