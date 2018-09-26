#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <numeric>
#include <iterator>
#include <random>
#include <cmath>
#include <fstream>
#include <time.h>
#include <functional>

using std::cout;
using std::endl;
using std::vector;
using std::string;

// Problem Type
std::ifstream ifs("./TSP1.txt");
const string problem_type = "MIN";          // "MIN" or "MAX"

// Read data
const int N_PTS = 1000;
vector<vector<double>> point_table(N_PTS, vector<double>(2));

// Results
int const MAX_EVAL=100000;
int num_eval;
int PATH_OUT_FREQ = MAX_EVAL;
std::ofstream ofs_history("./history.csv");
std::ofstream ofs_bestpath("./bestpath.csv");

// Random Seed
std::random_device seed_gen;
std::mt19937 mt_engine(seed_gen());


/* Create random path */
vector<int> create_initial_solution(){

    vector<int> vec(N_PTS);

    //Ordered Vector
    for (int i=0; i<N_PTS; ++i){vec[i]=i;};

    // Shuffle
    std::shuffle(vec.begin(), vec.end(), mt_engine);

    return vec;
}

/* Evaluation */
double fitness_func(vector<int> const &path){

    double distance = 0.0;

    for(int i=0; i<N_PTS; i++){

        int index0 = path[i];
        int index1 = path[i + 1];
        if(i==N_PTS-1) index1 = path[0];

        double x0 = point_table[index0][0];
        double y0 = point_table[index0][1];
        double x1 = point_table[index1][0];
        double y1 = point_table[index1][1];

        distance = distance +
                   std::sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
    }

    if (problem_type=="MIN"){
        return 1.0/distance;
    }else if(problem_type=="MAX"){
        return distance;
    }

    return 0.0;
}

/* Output function */
void output(vector<int> const &best_ind, double const &best_score){

    double best_path[N_PTS][2];

    ofs_history << num_eval << ", " << best_score << endl;

    if ( (num_eval + 1) % PATH_OUT_FREQ == 0){
        for (int i=0; i<N_PTS; i++){
            best_path[i][0] = point_table[best_ind[i]][0];
            best_path[i][1] = point_table[best_ind[i]][1];
            ofs_bestpath << best_path[i][0] << ", " << best_path[i][1] << endl;
        }
        ofs_bestpath << best_path[0][0] << ", " << best_path[0][1] << endl;
    }
}


int main(){

    // Read data.
    std::string str;
    int i=0;
    while(getline(ifs, str)) {
		sscanf(str.data(), "%lf,%lf", &point_table[i][0], &point_table[i][1]);
        i++;
        if (i==N_PTS) break;
    }

    vector<int> best_solution;
    vector<int> solution;
    double score;
    double best_score;
    double distance;

    // Initialize
    best_solution = create_initial_solution();
    best_score = fitness_func(best_solution);

    for (num_eval=0; num_eval < MAX_EVAL; num_eval++){

        // New solution
        solution = create_initial_solution();
        score = fitness_func(solution);

        // Update.
        if (score > best_score){
            best_solution = solution;
            best_score = score;
        }

        // --- Output
        if (problem_type=="MAX"){
            distance = best_score;
        }else if(problem_type=="MIN"){
            distance = 1.0/best_score;
        }

        output(best_solution, best_score);
        cout << "Generation: " << num_eval << " / Best Fitness: " << distance << endl;
    }

    return 0;
}
