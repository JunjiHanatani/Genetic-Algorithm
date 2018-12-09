#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iostream>
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
int N_TETRA;
int N_SYMMETRIC_PAIR;

vector<Mass> mass;
vector<Spring> spring;

vector<vector<double>> nodeInitialPosition;
vector<double> springInitialRestlength;
vector<double> breathe_amp;
vector<double> breathe_phase;
vector<double> breathe_offset;
vector<std::array<int, 8>> vertices;
vector<std::array<int, 4>> tetra_vertices;
vector<std::array<int, 6>> octa_vertices;
vector<vector<int>> symmetric_pair;
vector<double> range_min;
vector<double> range_max;
double robot_diameter;

// Node parameter
const double nodeRadius = 0.008;
const double nodeMass = 0.1;
vector<double> nodePositionOffset = {0.0, 0.0, 0.0};

// Edge parameter
double UnitLength = 0.1;
const double springK = 5e3;

// Actuator
double initial_breathe_amp = 0.0;
double initial_breathe_phase = 0.0;
double initial_breathe_offset = 0.0;
const double breathe_omega = 2.0*PI*1.0;

// Cube structure.
static const int N=4;
int arr[N][N][N];

//int arr[N][N][N] = {{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
//                    {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
//                    {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}};

//int arr[N][N][N] = {{{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}},
//                    {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}},
//                    {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}},
//                    {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}};

void GenerateCube(void){

  std::fill(arr[0][0], arr[N][0], 1);

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

  // --- mass position.
  for (int i=0; i<N_MASS; i++){
    nodeInitialPosition.push_back(mass[i].p);
  }

  range_min = {1e6, 1e6, 1e6};
  range_max = {-1e6, -1e6, -1e6};
  for (int i=0; i<N_MASS; i++){
    for (int j=0; j<3; j++){
      if(range_min[j] > mass[i].p[j]){
        range_min[j] = mass[i].p[j];
      }else if(range_max[j] < mass[i].p[j]){
        range_max[j] = mass[i].p[j];
      }
    }
  }
  //robot_diameter = calcDistance(range_min, range_max);
  robot_diameter = 0.01;
  range_min = {-0.22, 0.0,  0.0};
  range_max = { 0.22, 0.0,  0.4};

  FindSymmetricSpring();
  CubeLog();
}

void GenerateTetra(void){

  N_TETRA = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){

        if(arr[i][j][k]==1){
            N_TETRA += 1;

            /* ----- Generate Unit Tetrahedron ----- */

            double x=0.0;
            if ((j+k)%2==0) {
              x = (double) i;
            }else if((j+k)%2==1){
              x = (double) i + 0.5;
            }

            double y=0.0;
            if (k%2==0) {
              y = (double) j * sqrt(3.0)/2.0;
            }else if(k%2==1){
              y = (double) j * sqrt(3.0)/2.0 + sqrt(3.0)/6.0;
            }

            double z = k * sqrt(2.0/3.0);

            double pt[4][3] = {{x,     y,               z              },
                               {x+1.0, y,               z              },
                               {x+0.5, y+sqrt(3.0)/2.0, z              },
                               {x+0.5, y+sqrt(3.0)/6.0, z+sqrt(2.0/3.0)}};

            std::array<int, 4> vertex;

            for (int l=0; l<4; l++){

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

            tetra_vertices.push_back(vertex);

            /* ----- End of Unit Cube Generation----- */
        }
      }
    }
  }

  N_MASS = mass.size();

  for (int i=-1; i<N+1; i++){
    for (int j=-1; j<N+1; j++){
      for (int k=-1; k<N+1; k++){

        /* ----- Generate Unit Tetrahedron ----- */

        double x=0.0;
        if ((j+k)%2==0) {
          x = (double) i + 1.0;
        }else if((j+k)%2==1){
          x = (double) i + 0.5;
        }

        double y=0.0;
        if (k%2==0) {
          y = (double) j * sqrt(3.0)/2.0;
        }else if(k%2==1){
          y = (double) j * sqrt(3.0)/2.0 + sqrt(3.0)/6.0;
        }

        double z = k * sqrt(2.0/3.0);

        double pt[6][3] = {{x,     y                   , z              },
                           {x+0.5, y+sqrt(3.0)/6.0     , z+sqrt(2.0/3.0)},
                           {x+0.5, y+sqrt(3.0)/2.0     , z              },
                           {x    , y+sqrt(3.0)*2.0/3.0 , z+sqrt(2.0/3.0)},
                           {x-0.5, y+sqrt(3.0)/2.0     , z              },
                           {x-0.5, y+sqrt(3.0)/6.0     , z+sqrt(2.0/3.0)}};

        std::array<int, 6> vertex{-1, -1, -1, -1, -1, -1};

        for (int l=0; l<6; l++){

          for(int id=0; id<N_MASS; id++){
            if (mass[id].p[0] == pt[l][0] && mass[id].p[1] == pt[l][1] && mass[id].p[2] == pt[l][2]){
              vertex[l] = id;
              break;
            }
          }

        }

        octa_vertices.push_back(vertex);

      }
    }
  }

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
      if (length <= UnitLength*sqrt(2)+1e-2){
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

  // --- Store initial mass position.

  for (int i=0; i<N_MASS; i++){
    nodeInitialPosition.push_back(mass[i].p);
  }

  range_min = {1e6, 1e6, 1e6};
  range_max = {-1e6, -1e6, -1e6};
  for (int i=0; i<N_MASS; i++){
    for (int j=0; j<3; j++){
      if(range_min[j] > mass[i].p[j]){
        range_min[j] = mass[i].p[j];
      }else if(range_max[j] < mass[i].p[j]){
        range_max[j] = mass[i].p[j];
      }
    }
  }
  robot_diameter = calcDistance(range_min, range_max);

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
    mass[i].p = add(nodeInitialPosition[i], nodePositionOffset);
    mass[i].v = {0.0, 0.0, 0.0};
    mass[i].a = {0.0, 0.0, 0.0};
  }

  _Breathe = true;

  t = -1.0;
  nt = 0;

}

vector<int> Decorder(vector<vector<double>> parameter_table, string REPRESENTATION){

  int num = parameter_table.size();
  for (int i=0; i<N_SPRING; i++) spring[i].k = 0.0;
  for (int i=0; i<N_MASS; i++) mass[i].m = 0.0;
  int num_of_valid_masses = 0;
  int num_of_valid_springs = 0;

  /* DIRECT REPRESENTATION */
  if(REPRESENTATION == "direct"){

    for(int i=0; i<N_SPRING; i++){
      if (parameter_table[i][4]==1.0){
        breathe_offset[i] = parameter_table[i][0];
        breathe_amp[i] = parameter_table[i][1];
        breathe_phase[i] = parameter_table[i][2];
        spring[i].k = springK * pow(10, parameter_table[i][3]);
        mass[spring[i].masses[0]].m = nodeMass;
        mass[spring[i].masses[1]].m = nodeMass;
        num_of_valid_springs += 1;
      }
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
          if(parameter_table[i][4]==1.0){
            breathe_offset[j] += parameter_table[i][0];
            breathe_amp[j] += parameter_table[i][1];
            breathe_phase[j] += parameter_table[i][2];
            spring[i].k = springK * pow(10, parameter_table[i][3]);
            mass[spring[i].masses[0]].m = nodeMass;
            mass[spring[i].masses[1]].m = nodeMass;
            num_of_valid_springs += 1;
          }
        }
      }
    }

  }

  /* SYMMETRIC REPRESENTATION */
  else if (REPRESENTATION =="symmetric"){

    for (int i=0; i<N_SYMMETRIC_PAIR; i++){

      int id1 = symmetric_pair[i][0];
      int id2 = symmetric_pair[i][1];

      if(parameter_table[i+1][4]==1.0){
        breathe_offset[id1] = parameter_table[i+1][0];
        breathe_amp[id1] = parameter_table[i+1][1];
        breathe_phase[id1] = parameter_table[i+1][2];
        spring[id1].k = springK * pow(10, parameter_table[i+1][3]);
        mass[spring[id1].masses[0]].m = nodeMass;
        mass[spring[id1].masses[1]].m = nodeMass;

        breathe_offset[id2] = parameter_table[i+1][0];
        breathe_amp[id2] = parameter_table[i+1][1];
        breathe_phase[id2] = parameter_table[i+1][2] + parameter_table[0][2];
        spring[id2].k = springK * pow(10, parameter_table[i+1][3]);
        mass[spring[id2].masses[0]].m = nodeMass;
        mass[spring[id2].masses[1]].m = nodeMass;

        num_of_valid_springs += 2;
      }
    }

  }

  /* GENERATIVE REPRESENTATION */
  else if (REPRESENTATION =="generative"){

    for (int i=0; i<N_SPRING; i++){

      vector<int> sphere_list;
      breathe_offset[i] = 0.0;
      breathe_amp[i] = 0.0;
      breathe_phase[i] = 0.0;
      vector<double> k_factor(N_SPRING, 0.0);

      bool flag0=false; bool flag1=false;

      for (int j=0; j<num; j++){

        double x = parameter_table[j][4];
        double y = parameter_table[j][5];
        double z = parameter_table[j][6];
        double r = parameter_table[j][7];

        vector<double> pt = {x, y, z};
        vector<double> pt1 = mass[spring[i].masses[0]].p;
        vector<double> pt2 = mass[spring[i].masses[1]].p;

        if(calcDistance(pt, pt1) < r){
          sphere_list.push_back(j);
          flag0 = true;
        }

        if(calcDistance(pt, pt2) < r){
          sphere_list.push_back(j);
          flag1 = true;
        }
      }

      int list_size = sphere_list.size();

      if(flag0 && flag1){
        for (int id: sphere_list){
          breathe_offset[i] += parameter_table[id][0];
          breathe_amp[i] += parameter_table[id][1];
          breathe_phase[i] += parameter_table[id][2];
          k_factor[i] += parameter_table[i+1][3];
        }
        breathe_offset[i] = breathe_offset[i]/list_size;
        breathe_amp[i] = breathe_amp[i]/list_size;
        breathe_phase[i] = breathe_phase[i]/list_size;
        spring[i].k = springK * pow(10, k_factor[i]/list_size);
        mass[spring[i].masses[0]].m = nodeMass;
        mass[spring[i].masses[1]].m = nodeMass;
        num_of_valid_springs += 1;
      }

    }

  }

  /* GENERATIVE REPRESENTATION 2*/
  else if (REPRESENTATION =="generative2"){

    Eigen::Matrix3d sigma[num];
    Eigen::Vector3d mean[num];

    for(int id=0; id<num; id++){

      mean[id] << parameter_table[id][4],
                  parameter_table[id][5],
                  parameter_table[id][6];

      sigma[id] << parameter_table[id][7],  parameter_table[id][10], parameter_table[id][12],
                   parameter_table[id][10], parameter_table[id][8],  parameter_table[id][11],
                   parameter_table[id][12], parameter_table[id][11], parameter_table[id][9];

    }

    Mvn mvn[num];
    for (int i=0; i<num; i++){
      mvn[i].Set(mean[i], sigma[i]);
    }


    for (double thresh=0.0; thresh<50.0; thresh+=0.1){

      num_of_valid_springs = 0;
      vector<double> k_factor(N_SPRING, 0.0);
      for (int i=0; i<N_SPRING; i++) spring[i].k = 0.0;
      for (int i=0; i<N_MASS; i++) mass[i].m = 0.0;

      for (int i=0; i<N_SPRING; i++){

        breathe_offset[i] = 0.0;
        breathe_amp[i] = 0.0;
        breathe_phase[i] = 0.0;

        Mass m1 = mass[spring[i].masses[0]];
        Mass m2 = mass[spring[i].masses[1]];

        Eigen::Vector3d pt1;
        Eigen::Vector3d pt2;
        pt1 << m1.p[0], m1.p[1], m1.p[2];
        pt2 << m2.p[0], m2.p[1], m2.p[2];

        vector<double> weight1(num);
        vector<double> weight2(num);
        for (int id=0; id<num; id++) weight1[id] = mvn[id].pdf(pt1);
        for (int id=0; id<num; id++) weight2[id] = mvn[id].pdf(pt2);

        double sum1 = std::accumulate(weight1.begin(), weight1.end(), 0.0);
        double sum2 = std::accumulate(weight2.begin(), weight2.end(), 0.0);

        bool flag0=false;
        bool flag1=false;
        if(thresh < sum1 / (double)num) flag0 = true;
        if(thresh < sum2 / (double)num) flag1 = true;

        if(flag0 && flag1){
          double norm = sum1 + sum2;
          for (int id=0; id<num; id++){
            breathe_offset[i] += parameter_table[id][0] * weight1[id]/norm;
            breathe_offset[i] += parameter_table[id][0] * weight2[id]/norm;
            breathe_amp[i] += parameter_table[id][1] * weight1[id]/norm;
            breathe_amp[i] += parameter_table[id][1] * weight2[id]/norm;
            breathe_phase[i] += parameter_table[id][2] * weight1[id]/norm;
            breathe_phase[i] += parameter_table[id][2] * weight2[id]/norm;
            k_factor[i] += parameter_table[id][3] * weight1[id]/norm;
            k_factor[i] += parameter_table[id][3] * weight2[id]/norm;
          }

          spring[i].k = springK * pow(10, k_factor[i]);
          mass[spring[i].masses[0]].m = nodeMass;
          mass[spring[i].masses[1]].m = nodeMass;
          num_of_valid_springs += 1;
        }

      }

      if (num_of_valid_springs < 200) {
        break;
      }

    }
  }

  for (Mass m: mass) if(m.m!=0.0) num_of_valid_masses += 1;
  vector<int> num_of_parts = {num_of_valid_masses, num_of_valid_springs};
  return num_of_parts;

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
