#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <mpi.h>
#include "utility.h"
#include "GeneticOperators.h"
#include "Vector3D.h"
#include "CubeGenerator.h"
#include "RecordLog.h"
#include "PhysicsEngine.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;

int gen = 0;
int N_GEN = 10;
static const int N_POP = 10;
static const int FREQ_OUT = 1;
static const string REPRESENTATION="symmetric";

static const int N_ELITES = 1;
static const int N_TOURNAMENT = 3;
static const int N_LAYERS = 5;
static const double P_MUT = 0.7;
static const double RATE_OVERAGE = 0.2;
int MAX_AGE[N_LAYERS] = {3, 6, 12, 24, 48};
static double const T_MAX = 9.0;
static double const T_MIN = 2.0;
double offset_range[2] = {-0.03, 0.03};
double amp_range[2] = {0.0, 0.05};
double phase_range[2] = {-PI, PI};
static const double thresh = 0.001;
static const double alpha = 2;

double duration[5];
static std::ofstream ofs_log;
int SIZE_OF_CHROMOSOME;
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
        ind.para.push_back({offset, amp, phase});
      }
      pop[i] = ind;
    }
    return pop;
}


// -----------------------------------------------
// EVALUATION
// -----------------------------------------------

void evaluate(vector<Individual>&pop){
  int num = pop.size();

  for (int i=0; i<num; i++) {
    InitializeCube();

    // Set breathe parameter.
    SetBreathe(pop[i].para, REPRESENTATION);

    vector<double> center_z;
    for (int nt=0; nt<NT_MAX; nt++){
      PhysicsEngine();

      if (nt%2000==0 && nt>=6000) {
        double sum=0.0;
        for(Mass m:mass) sum += m.p[2];
        center_z.push_back(sum/N_MASS);
      }
      t += dt;
    }

    pop[i].distance = calcNorm(mass[0].p);

    double ave = std::accumulate(center_z.begin(), center_z.end(), 0.0) / center_z.size();
    double err=0.0;
    for(double x:center_z) {
      err += (x - ave) * (x - ave);
    }
    pop[i].err = sqrt(err);

    double factor;
    if (thresh > pop[i].err){
        factor = 1.0;
    }else{
        factor = std::pow(thresh/pop[i].err, alpha);
    }

    pop[i].fitness = pop[i].distance * factor;

  }
}

void evaluateMPI(vector<Individual>&pop){

  int size, rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  double local_fitness[N_POP*3];
  double global_fitness[N_POP*3];
  double local_distance[N_POP*3];
  double global_distance[N_POP*3];
  double local_err[N_POP*3];
  double global_err[N_POP*3];
  double local_trajectory[N_POP*3*3*SIZE_TRAJECTORY];
  double global_trajectory[N_POP*3*3*SIZE_TRAJECTORY];

  int num = pop.size();
  int num_mpi = std::ceil((double)num/(double)size);

  for(int iter=0; iter<num_mpi; iter++){

    int i = rank*num_mpi + iter;

    if (i < num){

      InitializeCube();

      // Set breathe parameter.
      SetBreathe(pop[i].para, REPRESENTATION);

      // Begin Physics Simulation.
      vector<vector<double>> center_list;
      for (nt=0; nt<=NT_MAX; nt++){
        PhysicsEngine();

        if (nt%NT_CYCLE==0) {
          vector<double> cp= {0.0, 0.0, 0.0};
          for(Mass m:mass) cp = add(cp, m.p);
          cp = scaling(cp, (double)1.0/N_MASS);
          center_list.push_back(cp);
        }

        t += dt;
      }

      //local_distance[iter] = calcDistance(center_list[0], center_list.back());
      local_distance[iter] = fabs(center_list[0][0] - center_list.back()[0]);

      /* --- Calc. Error --- */
      vector<double> z_list;
      for (int j=AVE_RANGE_MIN; j<SIZE_TRAJECTORY; j++){
        z_list.push_back(center_list[j][2]);
      }
      double ave = std::accumulate(z_list.begin(), z_list.end(), 0.0) / z_list.size();
      double err=0.0;
      for(double z:z_list) err += (z - ave) * (z - ave);
      local_err[iter] = sqrt(err);

      /* --- Calc. Factor --- */
      double factor;
      if (thresh > local_err[iter]){
          factor = 1.0;
      }else{
          factor = std::pow(thresh/local_err[iter], alpha);
      }

      /* --- Fitness --- */
      local_fitness[iter] = local_distance[iter] * factor;

      /* --- Trajectory --- */
      for (int j=0; j<SIZE_TRAJECTORY; j++){
        for (int k=0; k<3; k++){
          local_trajectory[iter*SIZE_TRAJECTORY*3 + j*3 + k] = center_list[j][k];
        }
      }

    }
  }

  MPI_Allgather( &local_err, num_mpi, MPI_DOUBLE, &global_err, num_mpi, MPI_DOUBLE, MPI_COMM_WORLD );
  MPI_Allgather( &local_distance, num_mpi, MPI_DOUBLE, &global_distance, num_mpi, MPI_DOUBLE, MPI_COMM_WORLD );
  MPI_Allgather( &local_fitness, num_mpi, MPI_DOUBLE, &global_fitness, num_mpi, MPI_DOUBLE, MPI_COMM_WORLD );
  MPI_Allgather( &local_trajectory, num_mpi*SIZE_TRAJECTORY*3, MPI_DOUBLE,
                 &global_trajectory, num_mpi*SIZE_TRAJECTORY*3, MPI_DOUBLE, MPI_COMM_WORLD );

  for (int i=0; i<num; i++){
    pop[i].fitness = global_fitness[i];
    pop[i].distance = global_distance[i];
    pop[i].err = global_err[i];
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
    vector<Individual> elites(n_elites);
    for (int i=0; i<N_ELITES; i++) elites[i] = pop[i];
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
    vector<Individual> elites(N_ELITES);

    for (int id=0; id<N_LAYERS; id++){

        if (id != N_LAYERS-1){
            overages[id+1] = overageSelection(pops[id], id);
        }

        int M = N_POP*(1.0 - RATE_OVERAGE);
        int N = N_POP*RATE_OVERAGE;
        int m = pops[id].size();
        int n = overages[id].size();

        if (m < M && n < N){
            pops[id].insert(pops[id].end(), overages[id].begin(), overages[id].end());

        }else if(m < M && n > N){
            N = N_POP - m;
            if (N < n){
                elites = elitistSelection(overages[id], N_ELITES);
                //overages[id] = rouletteSelection(overages[id], N - N_ELITES);
                overages[id] = tournamentSelection(overages[id], N - N_ELITES, N_TOURNAMENT);
                overages[id].insert(overages[id].end(), elites.begin(), elites.end());
            }
            pops[id].insert(pops[id].end(), overages[id].begin(), overages[id].end());

        }else if(m > M && n < N){
            M = N_POP - n;
            if (M < m){
                elites = elitistSelection(pops[id], N_ELITES);
                //pops[id] = rouletteSelection(pops[id], M - N_ELITES);
                pops[id] = tournamentSelection(pops[id], M - N_ELITES, N_TOURNAMENT);
                pops[id].insert(pops[id].end(), elites.begin(), elites.end());
            }
            pops[id].insert(pops[id].end(), overages[id].begin(), overages[id].end());

        }else if (m >= M && n >= N){
            elites = elitistSelection(pops[id], N_ELITES);
            //pops[id] = rouletteSelection(pops[id], M - N_ELITES);
            pops[id] = tournamentSelection(pops[id], M - N_ELITES, N_TOURNAMENT);
            pops[id].insert(pops[id].end(), elites.begin(), elites.end());

            elites = elitistSelection(overages[id], N_ELITES);
            //overages[id] = rouletteSelection(overages[id], N - N_ELITES);
            overages[id] = tournamentSelection(overages[id], N - N_ELITES, N_TOURNAMENT);
            overages[id].insert(overages[id].end(), elites.begin(), elites.end());

            pops[id].insert(pops[id].end(), overages[id].begin(), overages[id].end());
        }

        vector<double> age_list;
        int min_age; int max_age;
        if (id==0){min_age = 1;}else{min_age = MAX_AGE[id-1]+1;}
        max_age = MAX_AGE[id];
        for (Individual ind:pops[id]) age_list.push_back(ind.age);
        ofs_log << "    LAYER" << id << "|";
        for (int age=min_age; age<=max_age; age++){
            size_t n_count = std::count(age_list.begin(), age_list.end(), age);
            if (n_count!=0) ofs_log << " AGE " << age << "(" << n_count << ")";
        }
        ofs_log << endl;
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
    }
    ofs_log << "Done. (" << cnt << "/" << num << ")" << endl;
}

void mutNormal(Individual &ind){
    std::normal_distribution<> normal_dist_a(0.0, 0.01);
    std::normal_distribution<> normal_dist_b(0.0, 0.01);
    std::normal_distribution<> normal_dist_c(0.0, 0.1);
    int i = get_rand_range_int(0, SIZE_OF_CHROMOSOME-1);
    //for (int i=0; i<SIZE_OF_CHROMOSOME; i++){
      //if (get_rand_range_dbl(0.0, 1.0) < 0.3){
        ind.para[i][0] += normal_dist_a(mt_engine);
        ind.para[i][1] += normal_dist_b(mt_engine);
        ind.para[i][2] += normal_dist_c(mt_engine);
      //}
    //}
}

// -----------------------------------------------
// CROSSOVER
// -----------------------------------------------

void crossover(vector<Individual> &pop1, const vector<Individual> &pop2){

    int num = pop1.size();
    vector<Individual> add_pop;
    //std::shuffle(pop2.begin(), pop2.end(), mt_engine);
    for(int i=0; i<num; i++){
        vector<Individual> children = oneptcx(pop1[i], pop2[i]);
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


void alpsGA(vector<Individual>pops[]){
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    Timer tm;
    ofs_log << " ----- GENERATION: " << gen << " ----- " <<endl;

    for (int id=N_LAYERS-1; id>=0; id--){
        ofs_log << "  LAYER: " << id << endl;

        if(pops[id].size()==0){
            ofs_log << "    No Population: " << endl;
            continue;
        }

        vector<Individual> elites = elitistSelection(pops[id], N_ELITES);
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
        ofs_log << "best fitness=" << std::setprecision(6) << offspring[0].fitness << endl;

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
    if (pops[0].size()==0){
        pops[0] = createInitialPop(N_POP);
        evaluateMPI(pops[0]);
        ofs_log << "    New population created. " << endl;
    }

    // --- Sort.
    for (int i=0; i<N_LAYERS; i++) sort_pop(pops[i]);

    duration[3] += tm.elapsed();

    // ---------------------------------------
    // Output
    // ---------------------------------------
    if (rank==0){
      tm.restart();
      ofs_log << "  OUTPUT" << endl;
      RecordLog(pops);
      duration[4] += tm.elapsed(); tm.restart();
    }
    ofs_log << "END" << endl; ofs_log << endl;

}

// ---------------------------------------------------
// Genetic Algorithm
// ---------------------------------------------------
void EvolveCube(void)
{
    Timer tm;
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    vector<Individual> pops[N_LAYERS];

    if (REPRESENTATION=="direct"){
      SIZE_OF_CHROMOSOME = N_SPRING;
    }else if(REPRESENTATION=="cubic"){
      SIZE_OF_CHROMOSOME = N_CUBE;
    }else if(REPRESENTATION=="symmetric"){
      SIZE_OF_CHROMOSOME = N_SYMMETRIC_PAIR + 1;
    }

    // Create initial populations.
    if (rank==0) {
        ofs_log.open("./log/log.csv");
        ofs_log << " < START EVOLUTION >" <<endl;
        ofs_log <<endl;
        ofs_log << " ----- INITIAL POPULATION ----- " <<endl;
    }
    pops[0] = createInitialPop(N_POP);

    // Evaluation
    evaluateMPI(pops[0]);
    sort_pop(pops[0]);

    // Log
    if (rank==0){
      ofs_log << "  OUTPUT" << endl;
      RecordLog(pops);
      ofs_log << "END" << endl; ofs_log << endl;
    }

    // Start evolution.
    for (gen=1; gen<=N_GEN; gen++){
        alpsGA(pops);
    }

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

void Restart(void)
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    Timer tm;
    vector<Individual> pops[N_LAYERS];

    // Read restart data.
    vector<vector<double>>df;
    string filename = "./log/restart.csv";

    df = read_csv(filename, 7, -1);
    for (vector<double> data: df){
      int layer = data[0];
      int ID = data[1];
      int locus = data[2];

      if (locus==0) {
        Individual ind;
        ind.age = data[7];
        ind.fitness = data[8];
        pops[layer].push_back(ind);
      }

      vector<double> parameter ={data[3], data[4], data[5], data[6]};
      pops[layer][ID].para.push_back(parameter);
    }

    // Evaluation
    for (int i=0; i<N_LAYERS; i++){
      evaluateMPI(pops[i]);
      sort_pop(pops[i]);
    }

    df = read_csv(filename, 5, 5);
    int INI_GEN = (int)df[0][0];

    df = read_csv(filename, 1, 1);
    SIZE_OF_CHROMOSOME = df[0][6];

    // Log
    if (rank==0) {
      ofs_log.open("./log/log.csv", std::ios::app);
      ofs_log << " < RESTART >" << endl;
      ofs_log << "  N_POP: "    << df[0][0]    << "  |"
              << "  N_LAYERS: " << df[0][3] << "  |"
              << "  SIZE_OF_CHROMOSOME: " << df[0][6] << endl;
      for (int i=0; i<N_LAYERS; i++){
        ofs_log << "  LAYER" << i << "| ";
        ofs_log << "best fitness=" << std::setprecision(6) << pops[i][0].fitness << endl;
      }
      ofs_log << endl;
    }

    // Start evolution.
    for (gen=INI_GEN+1; gen<=N_GEN; gen++){
        alpsGA(pops);
    }

    // Log.
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

    vector<Individual> totalPOP;
    for (int i=0; i<N_LAYERS; i++){
        totalPOP.insert(totalPOP.end(), pops[i].begin(), pops[i].end());
    }
    sort_pop(totalPOP);

    // Size of the population.
    int num = totalPOP.size();

    // Output
    if (gen%FREQ_OUT==0){

        // Output File
        std::ofstream ofs_para("./log/para_" + std::to_string(gen) + ".csv");
        std::ofstream ofs_fitness("./log/fitness_" + std::to_string(gen) + ".csv");
        std::ofstream ofs_trajectory("./log/trajectory_" + std::to_string(gen) + ".csv");

        // Parameter
        ofs_para << std::setprecision(8)
                 << "offset" << "," << "amplitude" << "," << "phase" << endl;
        for (int i=0; i<num; i++){
            for (vector<double> vec: totalPOP[i].para){
                ofs_para << vec[0] << "," << vec[1] << "," << vec[2] << endl;
            }
        }

        // Fitness
        ofs_fitness << "fitness" << "," << "distance"<< "," << "err" << "," << "age" << endl;

        for (int i=0; i<num; i++){
            ofs_fitness << totalPOP[i].fitness << ","
                        << totalPOP[i].distance << ","
                        << totalPOP[i].err << ","
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
    }else{
      ofs_best.open("./log/best_fitness.csv", std::ios::app);
    }

    if (gen==0){
        ofs_best << "generation"<<",";
        for (int i=0;i<N_LAYERS;i++) ofs_best << "LAYER" + std::to_string(i) <<",";
        ofs_best << endl;
    }

    ofs_best << gen << ",";
    for (int i=0; i<N_LAYERS; i++){
        if (pops[i].size()!=0){
            ofs_best << pops[i][0].fitness << ",";
        }else{
            ofs_best << ",";
        }
    }
    ofs_best << endl;
    ofs_best.close();

    /* --- Restart --- */
    std::ofstream ofs_restart("./log/restart.csv");

    ofs_restart << "N_POP" << "," << "N_ELITES" << "," << "N_TOURNAMENT" << ","
                << "N_LAYERS" << "," << "P_MUT" << "," << "RATE_OVERAGE" << ","
                << "SIZE_OF_CHROMOSOME" << endl;
    ofs_restart << N_POP << "," << N_ELITES << "," << N_TOURNAMENT << ","
                << N_LAYERS << "," << P_MUT << "," << RATE_OVERAGE << ","
                << SIZE_OF_CHROMOSOME << endl;

    ofs_restart << "MAX_AGE" << endl;
    for(int i=0; i<N_LAYERS; i++) ofs_restart << MAX_AGE[i] << ",";
    ofs_restart << endl;

    ofs_restart << "Generation" << endl;
    ofs_restart << gen << endl;

    ofs_restart << "LAYER" << "," << "indID" << "," << "springID" << ","
                << "offset" << "," << "amp" << "," << "phase1" << ","
                << "age" <<  endl;

    for (int i=0; i<N_LAYERS; i++){
      int indID = 0;
      for(Individual ind: pops[i]){
        int springID=0;
        for(vector<double> vec: ind.para){
          ofs_restart << std::setprecision(16)
                      << i << "," << indID << "," << springID << ","
                      << vec[0] << "," << vec[1] << "," << vec[2] << ","
                      << ind.age << endl;
          springID+=1;
        }
        indID+=1;
      }
    }

    cout << "Generation: " << std::setw(3) << std::right << gen
    << "  Fitness: " << std::fixed << std::setprecision(6) << totalPOP[0].fitness
    << "  Distance: " << totalPOP[0].distance
    << "  Err: " << totalPOP[0].err
    << endl;

}

void setBestIndividual(void){
  Individual ind;
  int i = 0;
  int n = N_GEN;

  if (REPRESENTATION=="direct"){
    SIZE_OF_CHROMOSOME = N_SPRING;
  }else if(REPRESENTATION=="cubic"){
    SIZE_OF_CHROMOSOME = N_CUBE;
  }else if(REPRESENTATION=="symmetric"){
    SIZE_OF_CHROMOSOME = N_SYMMETRIC_PAIR + 1;
  }

  int INI = 1 + SIZE_OF_CHROMOSOME * i;
  int FIN = INI + SIZE_OF_CHROMOSOME - 1;
  string filename = "./log/para_" + std::to_string(n) + ".csv";
  vector<vector<double>> parameter_table = read_csv(filename, INI, FIN);
  SetBreathe(parameter_table, REPRESENTATION);

}
