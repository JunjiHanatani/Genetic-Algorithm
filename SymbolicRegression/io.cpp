#include <iostream>
#include <fstream>
#include <cmath>
#include "io.h"
#include "evaluation.h"
#include "genetic_operator.h"

//#include "matplotlibcpp.h"
//namespace plt = matplotlibcpp;

using std::cout;
using std::endl;
using std::vector;
using std::string;

// Input
int const DOWN_SAMPLE_STEP = 10;
int const N_PTS = 1000/DOWN_SAMPLE_STEP;
std::ifstream ifs("./sample/function1.csv");

// Output
int const DIST_OUT_FREQ = NGEN/10;
int const PATH_OUT_FREQ = NGEN;
std::ofstream ofs_history("./history.csv");
std::ofstream ofs_distribution("./distribution.csv");
std::ofstream ofs_bestpath("./bestpath.csv");

vector<vector<double>> point_table(N_PTS, vector<double>(2));

// Read data.
void read_data(){
    std::string str;
    int i=0; int j=0;
    while(getline(ifs, str)) {
		if (i%DOWN_SAMPLE_STEP==0){
            sscanf(str.data(), "%lf,%lf", &point_table[j][0], &point_table[j][1]);
            j++;
		}
        i++;
    }
}

// Output data.
void output(vector<Tree> pop, bool verbose){

    vector<double> x_vec(N_PTS);
    vector<double> y_vec_est(N_PTS);
    vector<double> y_vec_ref(N_PTS);

    Tree best_ind = pop[0];
    double best_fitness = best_ind.fitness;
    double worst_fitness = pop[N_POP-1].fitness;

    double range = worst_fitness - best_fitness;

    int max_size=0;
    for (int i=0; i<N_POP; i++){
        int heap_size = pop[i].nodes.size();
        if (max_size < heap_size){max_size = heap_size;}
    }


    for (int i=0; i<N_PTS; i++){
        x_vec[i] = point_table[i][0];
        y_vec_ref[i] = point_table[i][1];
    }
    y_vec_est = calcFunction(best_ind.nodes, x_vec);


    ofs_history << g << ", " << best_fitness << endl;

    if ( (g + 1) % PATH_OUT_FREQ == 0){
        for (int i=0; i<N_PTS; i++){
            ofs_bestpath << x_vec[i] << ", " << y_vec_est[i] << endl;
        }
    }

    if (g % DIST_OUT_FREQ == 0){
        for(int i=0; i<N_POP; i++){
            ofs_distribution << g << ", " << pop[i].fitness << endl;
        }
    }

    if (verbose){
        cout << "Generation: " << g
             << " / Fitness: " << best_fitness
             << " / Range: "   << range
             << " / max tree depth: " << log2(max_size)
             << endl;
    }


    string func_name = showFunction(best_ind, false);
    /*
    if ( (g + 1) % 10 == 0){
        plt::clf();
        plt::xlim(0, 10);
        plt::ylim(-1, 1);
        plt::plot(x_vec, y_vec_ref, "o");
        plt::named_plot(func_name, x_vec, y_vec_est);
        plt::legend();
        plt::pause(0.001);
        //plt::show();
        //getchar();
    }
    */
}

