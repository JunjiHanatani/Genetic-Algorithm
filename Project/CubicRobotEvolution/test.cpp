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

int size, rank;
MPI_Comm_rank( MPI_COMM_WORLD, &rank );
MPI_Comm_size( MPI_COMM_WORLD, &size );

int global_arr[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
int local_arr[10] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

MPI_Scatter(&global_arr, 3, MPI_INT, &local_arr, 3, MPI_INT, 0, MPI_COMM_WORLD);

cout  << rank << ": ";
for(int i=0; i<3; i++){
  cout << local_arr[i] << " ";
}
cout << endl;

}

