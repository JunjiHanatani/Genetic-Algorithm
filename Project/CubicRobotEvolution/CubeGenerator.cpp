#include <math.h>
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include "PhysicsEngine.h"
#include "CubeGenerator.h"
#include "OpenGL.h"
#include "Vector3D.h"
#include "utility.h"
using std::vector;
using std::cout;
using std::endl;

// Structure
int N_MASS;
int N_SPRING;
int N_CUBE;
int N_SYMMETRIC_PAIR;

vector<Mass> mass;
vector<Spring> spring;

vector<vector<double>> nodeInitialPosition;
vector<double> springInitialRestlength;
vector<double> breathe_amp;
vector<double> breathe_phase;
vector<double> breathe_offset;
vector<std::array<int, 8>> vertices;
vector<vector<int>> symmetric_pair;

// Node parameter
const double nodeRadius = 0.008;
const double nodeMass = 0.1;
vector<double> nodePositionOffset = {0.0, 0.0, 0.0};

// Edge parameter
double UnitLength = 0.1;
const double springK = 1e3;

// Actuator
double initial_breathe_amp = 0.0;
double initial_breathe_phase = 0.0;
double initial_breathe_offset = 0.0;
const double breathe_omega = 2.0*PI*1.0;

// Cube structure.
static const int N=3;
int arr[N][N][N] = {{{1, 1, 0}, {1, 1, 0}, {0, 0, 0}},
                    {{1, 1, 0}, {1, 1, 0}, {0, 0, 0}},
                    {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
//int arr[N][N][N] = {{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
//                    {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
//                    {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}};


void GenerateCube(void){

  N_CUBE = 0;
  for (double i=0; i<N; i++){
    for (double j=0; j<N; j++){
      for (double k=0; k<N; k++){

        if(arr[(int)i][(int)j][(int)k]==1){
            N_CUBE += 1;
            /* ----- Generate Unit Cube ----- */

            double pt[8][3] = {{i,   j,   k},
                               {i+1, j,   k},
                               {i+1, j+1, k},
                               {i  , j+1, k},
                               {i  , j,   k+1},
                               {i+1, j,   k+1},
                               {i+1, j+1, k+1},
                               {i  , j+1, k+1}};

            std::array<int, 8> vertex;

            for (int l=0; l<8; l++){

              // Check if the node is overlapped with the other nodes.
              bool overlap=false;
              int overlappedID = 0;
              if (i + j + k != 0){
                for(Mass m:mass){
                  if (m.p[0] == pt[l][0] && m.p[1] == pt[l][1] && m.p[2] == pt[l][2]){
                    overlap = true;
                    break;
                  }
                  overlappedID += 1;
                }
              }

              // If not overlapped, then add the new node to the list.
              if (overlap){
                vertex[l] = overlappedID;
                continue;
              }else{
                Mass m;
                m.p = {pt[l][0], pt[l][1], pt[l][2]};
                mass.push_back(m);
                vertex[l] = mass.size()-1;
              }

            }

            vertices.push_back(vertex);

            /* ----- End of Unit Cube Generation----- */
        }
      }
    }
  }

  N_MASS = mass.size();

  vector<double> centroid = {0.0, 0.0, 0.0};
  for(Mass m:mass) centroid = add(centroid, m.p);
  centroid = scaling(centroid, (double)1.0/N_MASS);

  // Set initial mass conditions.
  for (int i=0; i<N_MASS; i++){
      mass[i].m = nodeMass;
      mass[i].p[0] -= centroid[0];
      mass[i].p[1] -= centroid[1];
      mass[i].p = scaling(mass[i].p, UnitLength);
      mass[i].v = {0.0, 0.0, 0.0};
      mass[i].a = {0.0, 0.0, 0.0};
  }

  // Set initial spring conditions.
  for (int i=0; i<N_MASS; i++){
    for (int j=i+1; j<N_MASS; j++){
      double length = calcDistance(mass[i].p, mass[j].p);
      if (length <= UnitLength*sqrt(3)+1e-2){
        Spring s;
        s.masses = {i, j};
        s.l0 = length;
        s.k = springK;
        spring.push_back(s);
      }
    }
  }

  N_SPRING = spring.size();

  // Actuator

  /*  Store initial conditions  */
  // --- spring lengths.
  for (int i=0; i<N_SPRING; i++){
    springInitialRestlength.push_back(spring[i].l0);
    breathe_amp.push_back(initial_breathe_amp);
    breathe_phase.push_back(initial_breathe_phase);
    breathe_offset.push_back(initial_breathe_offset);
  }

  // --- spring lengths.
  for (int i=0; i<N_MASS; i++){
    nodeInitialPosition.push_back(mass[i].p);
  }

  FindSymmetricSpring();
  CubeLog();
}

void InitializeCube(void){

  // --- spring lengths.
  for (int i=0; i<N_SPRING; i++){
    spring[i].l0 = springInitialRestlength[i];
    breathe_offset[i] = 0.0;
    breathe_amp[i] = 0.0;
    breathe_phase[i] = 0.0;
  }

  // --- Initial Position.
  for (int i=0; i<N_MASS; i++){
    mass[i].p = nodeInitialPosition[i];
    mass[i].v = {0.0, 0.0, 0.0};
    mass[i].a = {0.0, 0.0, 0.0};
  }

  _Breathe = true;

  t = -1.0;
  nt = 0;

}

void SetBreathe(vector<vector<double>> parameter_table, string REPRESENTATION){

  int num = parameter_table.size();

  /* DIRECT REPRESENTATION */
  if(REPRESENTATION == "direct"){

    for(int i=0; i<N_SPRING; i++){
      breathe_offset[i] = parameter_table[i][0];
      breathe_amp[i] = parameter_table[i][1];
      breathe_phase[i] = parameter_table[i][2];
    }

  }

  /* CUBIC REPRESENTATION */
  else if(REPRESENTATION == "cubic"){

    for (int i=0; i<N_CUBE; i++){

      for (int j=0; j<N_SPRING; j++){

        bool flag0=false; bool flag1=false;

        for (int k=0; k<8; k++){
          int nodeID = vertices[i][k];
          if (spring[j].masses[0] == nodeID) flag0 = true;
          if (spring[j].masses[1] == nodeID) flag1 = true;
        }

        if (flag0 && flag1){
          breathe_offset[j] += parameter_table[i][0];
          breathe_amp[j] += parameter_table[i][1];
          breathe_phase[j] += parameter_table[i][2];
        }
      }
    }

  }

  /* SYMMETRIC REPRESENTATION */
  else if (REPRESENTATION =="symmetric"){

    for (int i=0; i<N_SYMMETRIC_PAIR; i++){

      int id1 = symmetric_pair[i][0];
      int id2 = symmetric_pair[i][1];

      breathe_offset[id1] = parameter_table[i+1][0];
      breathe_amp[id1] = parameter_table[i+1][1];
      breathe_phase[id1] = parameter_table[i+1][2];
      breathe_offset[id2] = parameter_table[i+1][0];
      breathe_amp[id2] = parameter_table[i+1][1];
      breathe_phase[id2] = parameter_table[i+1][2] + parameter_table[0][2];
    }

  }

   /* SYMMETRIC REPRESENTATION */
  else if (REPRESENTATION =="generative"){

    for (int i=0; i<N_SPRING; i++){

      vector<int> sphere_list;
      breathe_offset[i] = 0.0;
      breathe_amp[i] = 0.0;
      breathe_phase[i] = 0.0;

      for (int j=0; j<num; j++){

        double x = parameter_table[j][3];
        double y = parameter_table[j][4];
        double z = parameter_table[j][5];
        double r = parameter_table[j][6];

        vector<double> pt = {x, y, z};
        vector<double> pt1 = mass[spring[i].masses[0]].p;
        vector<double> pt2 = mass[spring[i].masses[1]].p;

        if(calcDistance(pt, pt1) < r){
          sphere_list.push_back(j);
        }

        if(calcDistance(pt, pt2) < r){
          sphere_list.push_back(j);
        }
      }

      int list_size = sphere_list.size();

      if(list_size!=0){
        for (int id: sphere_list){
          breathe_offset[i] += parameter_table[id][0];
          breathe_amp[i] += parameter_table[id][1];
          breathe_phase[i] += parameter_table[id][2];
        }
        breathe_offset[i] = breathe_offset[i]/list_size;
        breathe_amp[i] = breathe_amp[i]/list_size;
        breathe_phase[i] = breathe_phase[i]/list_size;
      }

    }

  }

}

void FindSymmetricSpring(void){

  double y = 0.0;
  N_SYMMETRIC_PAIR = 0;
  for (int i=0; i<N_SPRING; i++){
    vector<double> p1 = mass[spring[i].masses[0]].p;
    vector<double> p2 = mass[spring[i].masses[1]].p;

    for (int j=i; j<N_SPRING; j++){
      vector<double> q1 = mass[spring[j].masses[0]].p;
      vector<double> q2 = mass[spring[j].masses[1]].p;

      if((p1[0] == q1[0] && p1[1]-y == y-q1[1] && p1[2] == q1[2] &&
          p2[0] == q2[0] && p2[1]-y == y-q2[1] && p2[2] == q2[2]) ||
         (p1[0] == q2[0] && p1[1]-y == y-q2[1] && p1[2] == q2[2] &&
          p2[0] == q1[0] && p2[1]-y == y-q1[1] && p2[2] == q1[2]))

          {
          vector<int> pair = {i, j};
          symmetric_pair.push_back(pair);
          N_SYMMETRIC_PAIR += 1;
          }

    }
  }
}

void CubeLog(void){
  static std::ofstream ofs_cube("./log/cube.csv");

  ofs_cube << "Node parameter" << endl;
  ofs_cube << "nodeRadius = " << nodeRadius << endl;
  ofs_cube << "nodeMass = " << nodeMass << endl;
  //ofs_cube << "nodePositionOffset = ("
  //         << nodePositionOffset[0] << ","
  //         << nodePositionOffset[1] << ","
  //         << nodePositionOffset[2] << ")" << endl;

  ofs_cube << endl;

  ofs_cube << "Edge parameter" << endl;
  ofs_cube << "UnitLength = " << UnitLength << endl;
  ofs_cube << "springK = " << springK << endl;

  ofs_cube << endl;

  ofs_cube << "Actuation" << endl;
  ofs_cube << "initial_breathe_offset = " << initial_breathe_offset << endl;
  ofs_cube << "initial_breathe_amp = " << initial_breathe_amp << endl;
  ofs_cube << "initial_breathe_phase = " << initial_breathe_phase << endl;
  ofs_cube << "breathe_omega = " << breathe_omega << endl;

  ofs_cube << endl;

  ofs_cube << "Cube structure" << endl;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      ofs_cube <<  "(";
      for (int k=0; k<N; k++){
        ofs_cube << arr[i][j][k] << ", ";
      }
      ofs_cube <<  "), ";
    }
    ofs_cube << endl;
  }

  ofs_cube << endl;

  ofs_cube << "MASS" << endl;
  for (int i=0; i<N_MASS; i++){
    ofs_cube << i << " (" << mass[i].p[0] << "," << mass[i].p[1] << "," << mass[i].p[2] << ")" << endl;
  }

  ofs_cube << endl;

  ofs_cube << "SPRING" << endl;
  for (int i=0; i<N_SPRING; i++){
    ofs_cube << i << " (" << spring[i].masses[0] << "," << spring[i].masses[1] << "), " << spring[i].l0 <<endl;
  }
  ofs_cube << endl;

  ofs_cube << "Symmetric Pair" << endl;
  for (int i=0; i<N_SYMMETRIC_PAIR; i++){
    ofs_cube << i << " (" << symmetric_pair[i][0] << "," << symmetric_pair[i][1] << ") " <<endl;
  }
}
