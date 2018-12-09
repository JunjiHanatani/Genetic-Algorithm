#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <mpi.h>
#include <Eigen/Dense>
#include "utility.h"
#include "GeneticOperators.h"
#include "Vector3D.h"
#include "CubeGenerator.h"
#include "PhysicsEngine.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

int gen = 0;
int N_GEN = 500;
static const int N_POP = 30;
static const int FREQ_OUT = 1;

string REPRESENTATION;
int SIZE_OF_CHROMOSOME;
int SIZE_OF_GENE;

static const int N_AREA = 5;
static const int N_ELITES = 3;
static const int N_TOURNAMENT = 3;
static const int N_LAYERS = 6;
static const double P_MUT = 0.5;
static const double RATE_OVERAGE = 0.2;
int MAX_AGE[N_LAYERS-1] = {5, 10, 20, 40, 80};
static double const T_MAX = 9.0;
static double const T_MIN = 2.0;
double offset_range[2] = {-0.03, 0.03};
double amp_range[2] = {0.0, 0.05};
double phase_range[2] = {-PI, PI};
double kfactor_range[2] = {-1.0, 1.0};
static const double thresh = 0.01;
static const double alpha = 1;

double duration[5];
Timer tm_global;
static std::ofstream ofs_log;
static int const NT_MAX = (int) (T_MAX + 1.0)/dt;          // = 20000;
static int const NT_MIN = (int) (T_MIN + 1.0)/dt;          // = 6000;
static int const NT_CYCLE = (int)2.0*PI/breathe_omega/dt;  // = 2000
static int const SIZE_TRAJECTORY = NT_MAX / NT_CYCLE + 1;  // = 11;
static int const AVE_RANGE_MIN = NT_MIN / NT_CYCLE;        // = 3;

// -----------------------------------------------
// CREATE INITIAL POPULATION
// -----------------------------------------------

vector<Individual> createInitialPop(int num){

    vector<Individual> pop(num);
    //Create population
    for (int i=0; i<num; i++){

      Individual ind;
      for (int j=0; j<SIZE_OF_CHROMOSOME; j++){

          double offset = get_rand_range_dbl(offset_range[0], offset_range[1]);
          double amp = get_rand_range_dbl(amp_range[0], amp_range[1]);
          double phase = get_rand_range_dbl(phase_range[0], phase_range[1]);
          double kfactor = get_rand_range_dbl(kfactor_range[0], kfactor_range[1]);

        if(REPRESENTATION=="generative"){
          double x = get_rand_range_dbl(range_min[0], range_max[0]);
          double y = get_rand_range_dbl(range_min[1], range_max[1]);
          double z = get_rand_range_dbl(range_min[2], range_max[2]);
          double r = get_rand_range_dbl(0.0, robot_diameter/2.0);
          ind.para.push_back({offset, amp, phase, kfactor, x, y, z, r});

        }else if(REPRESENTATION=="generative2"){
          double x = get_rand_range_dbl(-0.2, 0.2);//(range_min[0], range_max[0]);
          double y = get_rand_range_dbl( 0.0, 0.0);//(range_min[1], range_max[1]);
          double z = get_rand_range_dbl( 0.0, 0.4);//(range_min[2], range_max[2]);
          // Eigen::Matrix3d m = Eigen::MatrixXd::Random(3,3) * sqrt(robot_diameter/6.0);
          // Eigen::Matrix3d sigma = m.transpose() * m;
          double sx = get_rand_range_dbl(0.0, 0.05); //sigma(0, 0);
          double sy = get_rand_range_dbl(0.0, 0.1); //sigma(1, 1);
          double sz = get_rand_range_dbl(0.0, 0.05); //sigma(2, 2);
          double sxy = 0.0;//sigma(0, 1);
          double syz = 0.0;//sigma(1, 2);
          double szx = 0.0;//sigma(0, 2);
          double thresh = 0.0;//get_rand_range_dbl(1.0, 3.0);
          ind.para.push_back({offset, amp, phase, kfactor, x, y, z, sx, sy, sz, sxy, syz, szx, thresh});
        }else{
          double exist = (double) get_rand_range_int(0, 1);
          ind.para.push_back({offset, amp, phase, kfactor, exist});
        }

      }
      pop[i] = ind;
    }

    return pop;
}


// -----------------------------------------------
// EVALUATION
// -----------------------------------------------

void evaluateMPI(vector<Individual>&pop){

  int size, rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  double local_distance[N_POP*3];
  double global_distance[N_POP*3];
  double local_err[N_POP*3];
  double global_err[N_POP*3];
  double local_diameter[N_POP*3];
  double global_diameter[N_POP*3];
  int local_num_of_masses[N_POP*3];
  int global_num_of_masses[N_POP*3];
  int local_num_of_springs[N_POP*3];
  int global_num_of_springs[N_POP*3];
  double local_trajectory [N_POP*3 * 3*SIZE_TRAJECTORY];
  double global_trajectory[N_POP*3 * 3*SIZE_TRAJECTORY];
  double local_parameter  [N_POP*3 * SIZE_OF_CHROMOSOME*SIZE_OF_GENE];
  double global_parameter [N_POP*3 * SIZE_OF_CHROMOSOME*SIZE_OF_GENE];

  int num;
  int num_mpi;

  /* --- Scatter parameter data in process 0 to others --- */

  if (rank==0){
    num = pop.size();
    num_mpi = std::ceil((double)num/(double)size);
    int cnt = 0;
    for (int i=0; i<num; i++){
      for (int j=0; j<SIZE_OF_CHROMOSOME; j++){
        for (int k=0; k<SIZE_OF_GENE; k++){
          global_parameter[cnt] = pop[i].para[j][k];
          cnt += 1;
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&num_mpi, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(&global_parameter, num_mpi*SIZE_OF_CHROMOSOME*SIZE_OF_GENE, MPI_DOUBLE,
              &local_parameter,  num_mpi*SIZE_OF_CHROMOSOME*SIZE_OF_GENE, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  /* --- Start parallel computing --- */
  for(int iter=0; iter<num_mpi; iter++){

    int i = rank*num_mpi + iter;

    if (i < num){
      /* --- Initialize the cube parameter --- */
      InitializeCube();

      /* --- Set breathe parameter. --- */
      vector<vector<double>> para_vec;
      for (int j=0; j<SIZE_OF_CHROMOSOME; j++){
        vector<double> gene;
        for (int k=0; k<SIZE_OF_GENE; k++){
          int id = SIZE_OF_GENE*SIZE_OF_CHROMOSOME*iter+j*SIZE_OF_GENE+k;
          gene.push_back(local_parameter[id]);
        }
        para_vec.push_back(gene);
      }

      /* --- Decording --- */
      vector<int> num_of_parts = Decorder(para_vec, REPRESENTATION);
      local_num_of_masses[iter] = num_of_parts[0];
      local_num_of_springs[iter] = num_of_parts[1];

      /* --- Calc. diameter --- */
      local_diameter[iter] = calcDiameter(mass);

      /* --- Begin physics simulation. --- */
      vector<vector<double>> center_list;
      for (nt=0; nt<=NT_MAX; nt++){
        PhysicsEngine();

        if (nt%NT_CYCLE==0) {
          vector<double> cp= {0.0, 0.0, 0.0};

          for(Mass m:mass) {
            if (m.m==0.0) continue;
            cp = add(cp, m.p);
          }
          if (num_of_parts[0]!=0) cp = scaling(cp, (double)1.0/num_of_parts[0]);
          center_list.push_back(cp);
        }

        t += dt;
      }

      /* --- Calculate Fitness --- */
      vector<double> results = calcTraveleLength(center_list);
      local_distance[iter] = results[0];
      local_err[iter] = results[1];

      /* --- Trajectory --- */
      for (int j=0; j<SIZE_TRAJECTORY; j++){
        for (int k=0; k<3; k++){
          local_trajectory[iter*SIZE_TRAJECTORY*3 + j*3 + k] = center_list[j][k];
        }
      }

    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather( &local_err, num_mpi, MPI_DOUBLE, &global_err, num_mpi, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Gather( &local_distance, num_mpi, MPI_DOUBLE, &global_distance, num_mpi, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Gather( &local_diameter, num_mpi, MPI_DOUBLE, &global_diameter, num_mpi, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Gather( &local_num_of_masses, num_mpi, MPI_INT, &global_num_of_masses, num_mpi, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Gather( &local_num_of_springs, num_mpi, MPI_INT, &global_num_of_springs, num_mpi, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Gather( &local_trajectory, num_mpi*SIZE_TRAJECTORY*3, MPI_DOUBLE,
                 &global_trajectory, num_mpi*SIZE_TRAJECTORY*3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  if(rank==0){
    for (int i=0; i<num; i++){
      pop[i].distance = global_distance[i];
      pop[i].diameter = global_diameter[i];
      pop[i].err = global_err[i];
      pop[i].num_of_masses = global_num_of_masses[i];
      pop[i].num_of_springs = global_num_of_springs[i];
      pop[i].fitness = calcFitness(pop[i].distance, pop[i].err, pop[i].diameter);
    }

    for (int i=0; i<num; i++){
      vector<vector<double>> trajectory;
      for (int j=0; j<SIZE_TRAJECTORY; j++){
        double x = global_trajectory[i*SIZE_TRAJECTORY*3 + j*3 + 0];
        double y = global_trajectory[i*SIZE_TRAJECTORY*3 + j*3 + 1];
        double z = global_trajectory[i*SIZE_TRAJECTORY*3 + j*3 + 2];
        trajectory.push_back({x, y, z});
      }
      pop[i].trajectory = trajectory;
    }
  }

  ofs_log << "Done." << endl;
}

double calcDiameter(vector<Mass>& mass){

  double diameter=0.0;
  int num_of_mass = 0;

  vector<double> range_min = {1e6, 1e6, 1e6};
  vector<double> range_max = {-1e6, -1e6, -1e6};

  for (int i=0; i<N_MASS; i++){
    if(mass[i].m==0.0) continue;
    for (int j=0; j<3; j++){
      if(range_min[j] > mass[i].p[j]){
        range_min[j] = mass[i].p[j];
      }else if(range_max[j] < mass[i].p[j]){
        range_max[j] = mass[i].p[j];
      }
    }
    num_of_mass += 1;
  }

  if (num_of_mass!=0) diameter = calcDistance(range_min, range_max);

  return diameter;
}

vector<double> calcTraveleLength(vector<vector<double>> &center_list){

  double distance=0.0;
  double err=0.0;
  vector<double> results;

  if (!std::isnan(center_list[0][0])){
    /* --- Calc. travel length --- */
    //distance = calcDistance(center_list[0], center_list.back());
    distance = center_list[SIZE_TRAJECTORY-1][0] - center_list[AVE_RANGE_MIN][0];

    /* --- Calc. error --- */
    vector<double> z_list, y_list, x_list;
    for (int j=AVE_RANGE_MIN; j<SIZE_TRAJECTORY; j++){
      x_list.push_back(center_list[j][0]);
      y_list.push_back(center_list[j][1]);
      z_list.push_back(center_list[j][2]);
    }

    int size_of_list = x_list.size();
    double z_ave = std::accumulate(z_list.begin(), z_list.end(), 0.0) / size_of_list;
    double y_ave = std::accumulate(y_list.begin(), y_list.end(), 0.0) / size_of_list;
    //double slope = y_list.back()/x_list.back();
    err = 0.0;
    for(int k=0; k<size_of_list; k++){
      //err += fabs(z_list[k] - ave);
      //err += sqrt(pow ((z_list[k] - z_ave), 2) + pow ((y_list[k] - x_list[k]*slope) * 0.2, 2));
      err += sqrt(pow ((z_list[k] - z_ave), 2) + pow ((y_list[k] - y_ave)* 0.2, 2));
    }
    err = err/size_of_list;

    /*
    vector<double> y_abs_list;
    for (int j=AVE_RANGE_MIN; j<SIZE_TRAJECTORY; j++){
      y_abs_list.push_back(fabs(center_list[j][1]));
    }
    err = *std::max_element(y_abs_list.begin(), y_abs_list.end());
    */
  }

  results = {distance, err};

  return results;
}

double calcFitness(double distance, double err, double diameter){

  double fitness=0.0;
  double factor=1.0;

  if (diameter!=0.0){
    /* Penalty */
    if (thresh < err){
        factor = std::pow(thresh/err, alpha);
    }

    /* Fitness */
    fitness = distance/diameter/T_MAX * 2.0*PI/breathe_omega*factor;
  }

  return fitness;
}

// -----------------------------------------------
//  SORT
// -----------------------------------------------

bool operator<(const Individual& left, const Individual& right){
  return left.fitness > right.fitness ;
}

void sort_pop(vector<Individual> &pop){
    std::sort(pop.begin(), pop.end());
}


// -----------------------------------------------
// SELECTION
// -----------------------------------------------

vector<Individual> tournamentSelection(vector<Individual> const &pop, int n_offspring, int n_tournament){

    vector<Individual> offspring(n_offspring);
    int rand_index;
    int min_index;

    for(int i=0; i<n_offspring; i++){
        min_index = pop.size();
        for (int j=0; j<n_tournament; j++){
            rand_index = get_rand_range_int(0, pop.size()-1);
            if (rand_index < min_index) min_index = rand_index;
        }
        offspring[i] = pop[min_index];
    }

    return offspring;
}

vector<Individual> rouletteSelection(vector<Individual> &pop, int n_offspring){

    vector<Individual> offspring;
    vector<double> rand_list(n_offspring);

    // Calculate sum of the fitness over all population.
    double sum_fitness = std::accumulate(pop.begin(), pop.end(), 0.0,
                     [](double sum, Individual& ind ){ return sum+ ind.fitness; } );

    // Generate random list.
    for (int i=0; i<n_offspring; i++){
        rand_list[i] = get_rand_range_dbl(0.0, sum_fitness);
    }

    // Sort random_list.
    std::sort(rand_list.begin(), rand_list.end());
    std::reverse(rand_list.begin(), rand_list.end());

    double thresh = 0.0;
    for (Individual ind: pop){
        thresh += ind.fitness;
        while(rand_list.size()!=0){
            double rand = rand_list.back();
            if (rand<thresh){
                offspring.push_back(ind);
                rand_list.pop_back();
            }else{
                break;
            }
        }
    }
    return offspring;
}

vector<Individual> elitistSelection(vector<Individual> const &pop, int n_elites){
    vector<Individual> elites;
    int num = pop.size();
    if (N_ELITES <= num){
      for (int i=0; i<N_ELITES; i++) elites.push_back(pop[i]);
    }else{
      elites = pop;
    }
    return elites;
}

vector<Individual> agebaseElitistSelection(vector<Individual> const &pop, int n_elites){

  vector<Individual> elites;

  for(int age=gen+1; age>0; age-=1){
    int cnt = 0;
    for(Individual ind: pop){
      if(ind.age==age) {
        cnt += 1;
        elites.push_back(ind);
        if (cnt==n_elites)break;
      }
    }
  }

  return elites;
}

vector<Individual> overageSelection(vector<Individual>&offspring, int id){
    vector<Individual> overages;
    int num=offspring.size();
    for (int i=num-1; i>=0; i--){
        if(offspring[i].age > MAX_AGE[id]){
            overages.insert(overages.begin(), offspring[i]);
            offspring.erase(offspring.begin()+i);
        }
    }
    return overages;
}

void agelayeredSelection(vector<Individual>pops[]){

    vector<Individual> overages[N_LAYERS];

    for (int id=0; id<N_LAYERS; id++){

        /* --- Overages --- */
        if (id != N_LAYERS-1){
            overages[id+1] = overageSelection(pops[id], id);
        }

        pops[id].insert(pops[id].end(), overages[id].begin(), overages[id].end());
        if (pops[id].size() <= N_POP) continue;

        /* --- Elite Selection --- */
        sort_pop(pops[id]);
        vector<Individual> elites = agebaseElitistSelection(pops[id], N_ELITES);

        /* --- Selection --- */

        int N = N_POP - elites.size();
        sort_pop(elites);
        if (N >= (int) N_POP*0.3){
          pops[id] = tournamentSelection(pops[id], N, N_TOURNAMENT);
        }else{
          N = (int) N_POP * 0.3;
          pops[id] = tournamentSelection(pops[id], N , N_TOURNAMENT);
          elites = tournamentSelection(elites, N_POP - N, N_TOURNAMENT);
        }
        pops[id].insert(pops[id].end(), elites.begin(), elites.end());

    }

}

// -----------------------------------------------
// MUTATION
// -----------------------------------------------

void mutation(vector<Individual> &pop){
    int num = pop.size();
    int cnt = 0;
    for (int i=0; i<num; i++){

        if (get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mutNormal(pop[i]);
            cnt += 1;
        }

        if (REPRESENTATION=="generative" || REPRESENTATION=="generative2"){
          if(get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mutSphere(pop[i]);
            cnt += 1;
          }
          if(get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mutRadius(pop[i]);
            cnt += 1;
          }
        }

        if (REPRESENTATION == "direct" || REPRESENTATION == "cubic" || REPRESENTATION == "symmetric"){
          if(get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mutExist(pop[i]);
            cnt += 1;
          }
        }
    }
    ofs_log << "Done. (" << cnt << "/" << num << ")" << endl;
}

void mutNormal(Individual &ind){
    std::normal_distribution<> normal_dist_1(0.0, 0.01);
    std::normal_distribution<> normal_dist_5(0.0, 0.05);
    std::normal_distribution<> normal_dist_10(0.0, 0.1);
    int i = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    ind.para[i][0] += normal_dist_1(mt_engine);
    ind.para[i][1] += normal_dist_1(mt_engine);
    ind.para[i][2] += normal_dist_5(mt_engine);
    ind.para[i][3] += normal_dist_10(mt_engine);
}

void mutSphere(Individual &ind){
    std::normal_distribution<> normal_dist(0.0, 0.03);
    int i = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    ind.para[i][4] += normal_dist(mt_engine);
    ind.para[i][5] += normal_dist(mt_engine);
    ind.para[i][6] += normal_dist(mt_engine);
}

void mutRadius(Individual &ind){
    int i = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);

    if(REPRESENTATION=="generative"){
      std::normal_distribution<> normal_dist(0.0, 0.03);
      ind.para[i][7] += normal_dist(mt_engine);

    }else if(REPRESENTATION=="generative2"){
      std::normal_distribution<> normal_dist(0.0, 0.002);
      ind.para[i][7] += normal_dist(mt_engine);
      ind.para[i][8] += normal_dist(mt_engine);
      ind.para[i][9] += normal_dist(mt_engine);
    }
}

void mutExist(Individual &ind){
    int i = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    ind.para[i][3]=0.0;
}

// -----------------------------------------------
// CROSSOVER
// -----------------------------------------------

void crossover(vector<Individual> &pop1, vector<Individual> &pop2){

    int num = pop1.size();
    vector<Individual> add_pop;
    std::shuffle(pop2.begin(), pop2.end(), mt_engine);
    for(int i=0; i<num; i++){
        vector<Individual> children = area_cx(pop1[i], pop2[i]);
        if (children.size()!=0){
            add_pop.insert(add_pop.end(), children.begin(), children.end());
        }
    }
    pop1 = add_pop;
    ofs_log << "Done. (" << num << "->" << add_pop.size() << ")" << endl;
}

// --- Single-point crossover
vector<Individual> oneptcx(Individual &ind1, const Individual &ind2){

    int cut_pt = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    vector<Individual> children;

    // Swap.
    vector<vector<double>> para1 = ind1.para;
    vector<vector<double>> para2 = ind2.para;
    vector<vector<double>> child_para1;
    vector<vector<double>> child_para2;

    child_para1.insert(child_para1.end(), para1.begin(), para1.begin()+cut_pt);
    child_para1.insert(child_para1.end(), para2.begin()+cut_pt, para2.end());
    child_para2.insert(child_para2.end(), para2.begin(), para2.begin()+cut_pt);
    child_para2.insert(child_para2.end(), para1.begin()+cut_pt, para1.end());
    int age = std::max(ind1.age, ind2.age);

    // New individuals.
    Individual child1; child1.para=child_para1; child1.age=age;
    Individual child2; child2.para=child_para2; child2.age=age;
    children = {child1, child2};

    //path_check(ind1, "cx0");
    //path_check(ind2, "cx1");
    //path_check(child1, "cx2");
    //path_check(child2, "cx3");

    return children;
}

vector<Individual> oneptswap(Individual &ind1, const Individual &ind2){

    int swap_pt1 = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    int swap_pt2 = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    vector<Individual> children;

    // Swap.
    vector<vector<double>> para1 = ind1.para;
    vector<vector<double>> para2 = ind2.para;
    vector<vector<double>> child_para1 = para1;
    vector<vector<double>> child_para2 = para2;

    child_para1[swap_pt1] = para2[swap_pt2];
    child_para2[swap_pt2] = para1[swap_pt1];
    int age = std::max(ind1.age, ind2.age);

    // New individuals.
    Individual child1; child1.para=child_para1; child1.age=age;
    Individual child2; child2.para=child_para2; child2.age=age;
    children = {child1, child2};

    return children;
}

vector<Individual> area_cx(Individual &ind1, const Individual &ind2){

    int rand_area = get_rand_range_int(0, N_AREA-1);
    int num = get_rand_range_int(1, N_AREA-2);

    vector<Individual> children;
    vector<double> distance_list1, distance_list2;
    vector<double> pt0 = {ind1.para[rand_area][3], ind1.para[rand_area][4], ind1.para[rand_area][5]};

    for (int j=0; j<N_AREA; j++){

      double x = ind1.para[j][3];
      double y = ind1.para[j][4];
      double z = ind1.para[j][5];
      vector<double> pt1 = {x, y, z};
      distance_list1.push_back(calcDistance(pt0, pt1));

    }

    for (int j=0; j<N_AREA; j++){

      double x = ind2.para[j][3];
      double y = ind2.para[j][4];
      double z = ind2.para[j][5];
      vector<double> pt1 = {x, y, z};
      distance_list2.push_back(calcDistance(pt0, pt1));

    }

    vector<int> indices1 = argsort(distance_list1);
    vector<int> indices2 = argsort(distance_list2);

    // Swap.
    vector<vector<double>> para1 = ind1.para;
    vector<vector<double>> para2 = ind2.para;
    vector<vector<double>> child_para1 = ind1.para ;
    vector<vector<double>> child_para2 = ind2.para;

    for (int i=0; i<num; i++){
      int index1 = indices1[i];
      int index2 = indices2[i];
      child_para1[index1] = para2[index2];
      child_para2[index2] = para1[index1];
    }
    int age = std::max(ind1.age, ind2.age);

    // New individuals.
    Individual child1; child1.para=child_para1; child1.age=age;
    Individual child2; child2.para=child_para2; child2.age=age;
    children = {child1, child2};

    return children;
}


void alpsGA(vector<Individual>pops[]){
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    Timer tm;
    ofs_log << " ----- GENERATION: " << gen << " ----- " <<endl;

    for (int id=N_LAYERS-1; id>=0; id--){

        ofs_log << "  LAYER: " << id << endl;

        int num = pops[id].size();
        MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(num==0){
            ofs_log << "    No Population: " << endl;
            continue;
        }

        //vector<Individual> elites = elitistSelection(pops[id], N_ELITES);
        vector<Individual> elites = agebaseElitistSelection(pops[id], N_ELITES);
        vector<Individual> offspring = pops[id];

        // --------------------------------------
        // Crossover
        // --------------------------------------
        tm.restart();
        ofs_log << "    CROSSOVER:  ";

        vector<Individual> partners;
        for (int i=0; i<=id; i++){
            partners.insert(partners.begin(), pops[i].begin(), pops[i].end());
        }

        crossover(offspring, partners);
        duration[0] += tm.elapsed();

        // --------------------------------------
        // Mutation
        // --------------------------------------
        tm.restart();
        ofs_log << "    MUTATION:   ";
        mutation(offspring);
        duration[1] += tm.elapsed();

        // --------------------------------------
        //  Evaluation
        // --------------------------------------
        tm.restart();
        ofs_log << "    EVALUATION: ";
        evaluateMPI(offspring);
        offspring.insert(offspring.end(), elites.begin(), elites.end());
        sort_pop(offspring);
        for(Individual &ind:offspring) ind.age += 1;
        pops[id] = offspring;

        duration[2] += tm.elapsed();

    }

    // --------------------------------------
    // Selection
    // --------------------------------------

    tm.restart();
    ofs_log << "  SELECTION" << endl;

    // --- Age-layered Selection
    agelayeredSelection(pops);

    // --- Create new population
    int num = pops[0].size();
    MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (num==0){
      if(rank==0) pops[0] = createInitialPop(N_POP);
      ofs_log << "    New population created. Evaluation ";
      evaluateMPI(pops[0]);
    }

    // --- Sort.
    for (int i=0; i<N_LAYERS; i++) sort_pop(pops[i]);

    /* --- Output next generations --- */
    for (int id=0; id<N_LAYERS; id++){
      ofs_log << "    LAYER" << id << "|";
      if(pops[id].size()!=0) ofs_log << " Best fitness=" << pops[id][0].fitness;
      vector<double> age_list;
      for (Individual ind:pops[id]) age_list.push_back(ind.age);
      for (int age=1; age<=gen+1; age++){
          size_t n_count = std::count(age_list.begin(), age_list.end(), age);
          if (n_count!=0) ofs_log << " AGE " << age << "(" << n_count << ") ";
      }
      ofs_log << endl;
    }

    duration[3] += tm.elapsed();

    // ---------------------------------------
    // Output
    // ---------------------------------------

    tm.restart();
    ofs_log << "  OUTPUT";
    if (rank==0) RecordLog(pops);
    duration[4] += tm.elapsed(); tm.restart();

}

// ---------------------------------------------------
// Genetic Algorithm
// ---------------------------------------------------
void EvolveCube(void){

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if(rank==0){
      ofs_log.open("./log/log.csv");
      ofs_log << " < START EVOLUTION >" <<endl;
    }

    /* ----------- Set Representation. ------------*/
    char rep;
    if (rank==0){
      cout << endl << "Please select the representation mode" << endl;
      cout << "(d) direct / (c) cubic / (s) symmetric / (g) generative2" << endl;
      std::cin >> rep;
      cout << endl;
    }
    MPI_Bcast(&rep, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    setRepresentation(rep);
    ofs_log << "  Representation: " << REPRESENTATION << endl <<endl;

    /* -------- Create initial population. --------*/

    ofs_log << " ----- INITIAL POPULATION ----- " << endl;
    vector<Individual> pops[N_LAYERS];
    if(rank==0) pops[0] = createInitialPop(N_POP);
    ofs_log << "Done." << endl;

    /* ----------------- Evaluation. ------------- */
    ofs_log << "  EVALUATION ";
    evaluateMPI(pops[0]);
    sort_pop(pops[0]);

    /* ----------------- Output ------------------ */
    ofs_log << "  OUTPUT";
    if (rank==0){
      cout << "REPRESENTATION: " << REPRESENTATION << endl;
      cout << "N_POP = " << N_POP << endl;
      cout << "N_LAYERS = " << N_LAYERS << endl;
      cout << "N_GEN = " << N_GEN << endl;
      cout << "SIZE_OF_CHROMOSOME = " << SIZE_OF_CHROMOSOME << endl;
      cout << "SIZE_OF_GENE = " << SIZE_OF_GENE << endl << endl;
      RecordLog(pops);
    }

    /* ----------------- Evolution ----------------*/
    for (gen=1; gen<=N_GEN; gen++){
        alpsGA(pops);
    }

    /* ----------------- Output ------------------ */
    if (rank==0){
      cout << endl;
      cout << " --- Average elapsed time --- " << endl;
      cout << "Crossover : " << duration[0]/N_GEN << endl;
      cout << "Mutation  : " << duration[1]/N_GEN << endl;
      cout << "Evaluation: " << duration[2]/N_GEN << endl;
      cout << "Selection : " << duration[3]/N_GEN << endl;
      cout << "Output    : " << duration[4]/N_GEN << endl;
      cout << endl;
    }

}

void Restart(void){

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    vector<Individual> pops[N_LAYERS];

    if(rank==0){
      ofs_log.open("./log/log.csv", std::ios::app);
      ofs_log << " < RESTART >" << endl;
    }

    /* Read restart data */
    string filename = "./log/restart.csv";
    int ini_gen = readRestartData(pops, filename);

    /* Evaluation */
    for (int i=0; i<N_LAYERS; i++){
      ofs_log << "  Evaluation (Layer" + std::to_string(i) + ")";
      evaluateMPI(pops[i]);
      sort_pop(pops[i]);
    }

    // Log
    if (rank==0) {
      for (int i=0; i<N_LAYERS; i++){
        ofs_log << "  LAYER" << i << "| ";
        if(pops[i].size()==0){
          ofs_log << "No Population: " << endl;
        }else{
          ofs_log << "best fitness=" << std::setprecision(6) << pops[i][0].fitness << endl;
        }
      }
      ofs_log << endl;
    }

    if (rank==0){
      cout << "REPRESENTATION: " << REPRESENTATION << endl;
      cout << "N_POP = " << N_POP << endl;
      cout << "N_LAYERS = " << N_LAYERS << endl;
      cout << "N_GEN = " << N_GEN << endl;
      cout << "SIZE_OF_CHROMOSOME = " << SIZE_OF_CHROMOSOME << endl;
      cout << "SIZE_OF_GENE = " << SIZE_OF_GENE << endl << endl;
    }

    /* Start evolution. */
    for (gen=ini_gen+1; gen<=N_GEN; gen++){
        alpsGA(pops);
    }

    /* Log. */
    if (rank==0){
      cout << endl;
      cout << " --- Average elapsed time --- " << endl;
      cout << "Crossover : " << duration[0]/N_GEN << endl;
      cout << "Mutation  : " << duration[1]/N_GEN << endl;
      cout << "Evaluation: " << duration[2]/N_GEN << endl;
      cout << "Selection : " << duration[3]/N_GEN << endl;
      cout << "Output    : " << duration[4]/N_GEN << endl;
      cout << endl;
    }

}

void RecordLog(vector<Individual>pops[]){
    int num; /* Size of population */

    /* --- Merge all layers --- */
    vector<Individual> totalPOP;
    for (int i=0; i<N_LAYERS; i++){
        totalPOP.insert(totalPOP.end(), pops[i].begin(), pops[i].end());
    }
    num = totalPOP.size();
    sort_pop(totalPOP);

    // Output
    if (gen%FREQ_OUT==0 || gen==N_GEN){

        // Output File
        std::ofstream ofs_para("./log/parameter/para_" + std::to_string(gen) + ".csv");
        std::ofstream ofs_fitness("./log/fitness/fitness_" + std::to_string(gen) + ".csv");
        std::ofstream ofs_trajectory("./log/trajectory/trajectory_" + std::to_string(gen) + ".csv");

        // Parameter
        ofs_para << REPRESENTATION << endl;
        ofs_para << "offset" << "," << "amplitude" << "," << "phase" << endl;
        ofs_para << std::setprecision(16);
        for (int i=0; i<num; i++){
            for (int j=0; j<SIZE_OF_CHROMOSOME; j++){
                for (int k=0; k<SIZE_OF_GENE; k++){
                    ofs_para << totalPOP[i].para[j][k] << ",";
                }
                ofs_para << endl;
            }
        }

        // Fitness
        ofs_fitness << "fitness" << "," <<"distance"<< "," << "err" << "," << "diameter" << ","
                    << "Num of masses" << "," << "Num of springs"<< "," << "age" << endl;

        for (int i=0; i<num; i++){
            ofs_fitness << totalPOP[i].fitness << ","
                        << totalPOP[i].distance << ","
                        << totalPOP[i].err << ","
                        << totalPOP[i].diameter << ","
                        << totalPOP[i].num_of_masses << ","
                        << totalPOP[i].num_of_springs << ","
                        << totalPOP[i].age << endl;
        }

        // Trajectory
        ofs_trajectory << "time"<< "," << "x" << "," << "y" << "," << "z" << endl;

        for (int i=0; i<num; i++){
          for (int j=0; j<SIZE_TRAJECTORY; j++){
            ofs_trajectory << (double) j - 1.0 << ","
                           << totalPOP[i].trajectory[j][0] << ","
                           << totalPOP[i].trajectory[j][1] << ","
                           << totalPOP[i].trajectory[j][2] << endl;
          }
        }
    }

    /* --- Best fitness --- */

    static std::ofstream ofs_best;
    if (gen==0){
      ofs_best.open("./log/best_fitness.csv");
      ofs_best << "generation"<<",";
      for (int i=0;i<N_LAYERS;i++) ofs_best << "LAYER" + std::to_string(i) <<",";
      ofs_best << endl;

    }else{
      ofs_best.open("./log/best_fitness.csv", std::ios::app);
      ofs_best << gen << ",";
      for (int i=0; i<N_LAYERS; i++){
          if (pops[i].size()!=0){
              ofs_best << pops[i][0].fitness << ",";
          }else{
              ofs_best << ",";
          }
      }
      ofs_best << endl;
    }
    ofs_best.close();

    /* --- Restart --- */
    std::ofstream ofs_restart("./log/restart.csv");

    ofs_restart << REPRESENTATION << endl;
    ofs_restart << "N_POP" << "," << "N_ELITES" << "," << "N_TOURNAMENT" << ","
                << "N_LAYERS" << "," << "P_MUT" << "," << "RATE_OVERAGE" << ","
                << "SIZE_OF_CHROMOSOME" << endl;
    ofs_restart << N_POP << "," << N_ELITES << "," << N_TOURNAMENT << ","
                << N_LAYERS << "," << P_MUT << "," << RATE_OVERAGE << ","
                << SIZE_OF_CHROMOSOME << endl;

    ofs_restart << "MAX_AGE" << endl;
    for(int age: MAX_AGE) ofs_restart << age << ",";
    ofs_restart << endl;

    ofs_restart << "Generation" << endl;
    ofs_restart << gen << endl;

    ofs_restart << "LAYER" << "," << "indID" << "," << "springID" << "," << "age" << ","
                << "offset" << "," << "amp" << "," << "phase" <<  endl;

    for (int i=0; i<N_LAYERS; i++){
      int indID = 0;
      for(Individual ind: pops[i]){
        int springID=0;
        for(vector<double> vec: ind.para){
          ofs_restart << std::setprecision(16)
                      << i << "," << indID << "," << springID << "," << ind.age;
          for (double x: vec){
            ofs_restart << "," << x;
          }
          ofs_restart << endl;
          springID+=1;
        }
        indID+=1;
      }
    }

    cout << "Generation: " << std::setw(3) << std::right << gen << "  |"
         << "  Fitness: " << std::fixed << std::setprecision(6) << totalPOP[0].fitness
         << "  Distance: " << totalPOP[0].distance
         << "  Err: " << totalPOP[0].err
         << "  Diameter: " << totalPOP[0].diameter << "  |"
         << "  Time: " << std::setw(6) << std::setprecision(3) << std::right << tm_global.elapsed() << " sec"
         << endl;
    tm_global.restart();

    ofs_log << "        Done." << endl << endl;
}

void setRepresentation(char rep){

  if ( rep=='d'){
    REPRESENTATION="direct";
    SIZE_OF_CHROMOSOME = N_SPRING;
    SIZE_OF_GENE = 5;

  }else if(rep=='c'){
    REPRESENTATION="cubic";
    SIZE_OF_CHROMOSOME = N_CUBE;
    SIZE_OF_GENE = 5;

  }else if(rep=='s'){
    REPRESENTATION="symmetric";
    SIZE_OF_CHROMOSOME = N_SYMMETRIC_PAIR + 1;
    SIZE_OF_GENE = 5;

  }else if (rep=='g'){
    REPRESENTATION="generative2";
    SIZE_OF_CHROMOSOME = N_AREA;
    SIZE_OF_GENE = 14;
  }
  //}else if (rep=='g'){
  //  REPRESENTATION="generative";
  //  SIZE_OF_CHROMOSOME = N_AREA;
  //  SIZE_OF_GENE = 7;
  //}
}

int readRestartData(vector<Individual> pops[], string filename){

    // Read restart data.
    vector<vector<double>>df;
    int INI_GEN;

    /* Representation */
    vector<vector<string>> rep = read_csv_string(filename, 0, 0);
    string str = rep[0][0];
    setRepresentation(str.c_str()[0]);

    /* Parameter data */
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if (rank==0){
      df = read_csv(filename, 8, -1);
      for (vector<double> data: df){
        int layer = data[0];
        int ID = data[1];
        int locus = data[2];

        if (locus==0) {
          Individual ind;
          ind.age = data[3];
          pops[layer].push_back(ind);
        }

        vector<double> parameter;
        parameter.insert(parameter.end(), data.begin()+4, data.end());
        pops[layer][ID].para.push_back(parameter);
      }
    }

    df = read_csv(filename, 6, 6);
    INI_GEN = (int)df[0][0];

    return INI_GEN;
}
