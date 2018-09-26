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
#include <string>

using std::cout;
using std::endl;
using std::vector;
using std::string;

// Problem Type
std::ifstream ifs("./TSP1.txt");
const string problem_type = "MIN";    // "MIN" or "MAX"

// Read data
const int N_PTS = 1000;
vector<vector<double>> point_table(N_PTS, vector<double>(2));

// Output
int MAX_EVALS = 10000000;
int PATH_OUT_FREQ = MAX_EVALS;
int num_evals = 0;
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

void output(vector<int> const &best_ind, double const &global_best, bool file_output){
    double best_path[N_PTS][2];

    if (file_output){

        ofs_history << num_evals << ", " << global_best << endl;

        if ( (num_evals + 1) % PATH_OUT_FREQ == 0){
            for (int i=0; i<N_PTS; i++){
                best_path[i][0] = point_table[best_ind[i]][0];
                best_path[i][1] = point_table[best_ind[i]][1];
                ofs_bestpath << best_path[i][0] << ", " << best_path[i][1] << endl;
            }
            ofs_bestpath << best_path[0][0] << ", " << best_path[0][1] << endl;
        }

    }

}

/* One step hill-climb */

void hillclimb(vector<int> &best_vec, double &best_score, bool &climbing_is_done){

    bool move_made = false;
    vector<int> vec = best_vec;

    for (int i=0; i<N_PTS; i++){
        for (int j=0; j<N_PTS; j++){

            vec = best_vec;
            // Exchange two cities.
            int tmp = vec[i];
            vec[i] = vec[j];
            vec[j] = tmp;

            // Evaluate the new vector.
            double score = fitness_func(vec);
            num_evals += 1;

            // if the best vector updates, move on next climb.
            if (score > best_score){
                best_vec = vec;
                best_score = score;
                move_made = true;
                break;
            }

            // if num_evals reaches the max number, the calculation ends.
            if (num_evals >= MAX_EVALS){
                climbing_is_done = true;
                break;
            }

        }
        if (move_made == true) break;
        if (climbing_is_done == true) break;
    }

    if (move_made == false){
        climbing_is_done = true;
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

    // Random restart hill-climbing
    vector<int> solution;
    double global_best = 0.0;
    double score; double distance;
    bool climbing_is_done;

    for (int i=0; i<100; i++){

        // Random restart
        solution = create_initial_solution();
        score = fitness_func(solution);
        climbing_is_done = false;

        // Hill-climbing start
        while (climbing_is_done==false){

            // One step climbing
            hillclimb(solution, score, climbing_is_done);

            // Output
            if (problem_type=="MAX"){
                distance = score;
            }else if(problem_type=="MIN"){
                distance = 1.0/score;
            }

            cout << "Evaluations: " << num_evals
                 << " / Best Fitness: " << distance
                 << " / Move: " << climbing_is_done
                 << endl;

            if(global_best < score){
                global_best = score;
                output(solution, global_best, true);
            }
        }

        if (num_evals >= MAX_EVALS) break;
    }

    return 0;
}

