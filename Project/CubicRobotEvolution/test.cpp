#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <vector>
#include <math.h>
#include "Vector3D.h"
#include "CubeGenerator.h"
#include "utility.h"
#include "test.h"
#include "GeneticOperators.h"
using std::cout;
using std::endl;
using std::vector;

void test(void){


Individual ind;
ind.para =
{{-0.017007, 0.00676091, -1.10495},
 {-0.0166607, 0.0193244, 2.5296},
 {0.00705909, 0.042877, -2.41781},
 {0.0091884, 0.00854548, -0.891256}};
gen = 1;
vector<Individual> pop= {ind};
evaluateMPI(pop);

}

