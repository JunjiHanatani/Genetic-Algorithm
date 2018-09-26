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

// Data
std::ifstream ifs("../sample/TSP1.txt");
const int N_PTS = 1000;
vector<vector<double>> point_table(N_PTS, vector<double>(2));

// GA type
const string problem_type = "MIN";          // "MIN" or "MAX"
const string selection_type = "tournament"; // "roulette" or "tournament"
const string variation_type = "improved";   // "simple" or "improved"

// GA parameter
const int NGEN = 1000;
const int N_POP = 100;
const int N_ELITES = 5;
const double P_MUT = 0.3;
const int N_TOURNAMENT = 3;
int g;

// Results
int DIST_OUT_FREQ = NGEN/100;
int PATH_OUT_FREQ = NGEN;
std::ofstream ofs_history("./history.csv");
std::ofstream ofs_distribution("./distribution.csv");
std::ofstream ofs_bestpath("./bestpath.csv");

vector<int> best_ind;
double best_fitness;
double best_path[N_PTS][2];
vector<double> history;
vector<double> distribution;

// Random Seed
std::random_device seed_gen;
std::mt19937 mt_engine(seed_gen());


int get_rand_range_int(int min_val, int max_val) {
    std::uniform_int_distribution<int> gen_rand_uni_int( min_val, max_val );
    return gen_rand_uni_int(mt_engine);
}

double get_rand_range_dbl(double min_val, double max_val) {
    std::uniform_real_distribution<double> gen_rand_uni_real( min_val, max_val );
    return gen_rand_uni_real(mt_engine);
}

int findIndex( vector<int> vec, int value ){
    vector<int>::iterator iter = std::find( vec.begin(), vec.end(), value);
    size_t index = std::distance( vec.begin(), iter );
    if(index == vec.size())
        {
            return -1;
        }
    return index;
}

/* Create random path */

vector<vector<int>> create_initial_pop(){

    vector<vector<int>> pop(N_POP, vector<int>(N_PTS));
    vector<int> vec(N_PTS);

    //Original Vector
    for (int i=0; i<N_PTS; ++i){vec[i]=i;};

    //Create population
    for (int i=0; i<N_POP; i++){
        // Shuffle
        // obtain a device-based seed:
        std::shuffle(vec.begin(), vec.end(), mt_engine);
        pop[i] = vec;
    }

    return pop;
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

void sort_pop(vector<vector<int>> &pop, vector<double> const &scores){

    int num = pop.size();
    vector<vector<int>> sorted_pop(num, vector<int>(N_PTS));
    vector<int> indices(num);

    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(),
              indices.end(),
              [&](int i, int j){return scores[i] > scores[j];}
    );

    for (int i=0; i<num; i++) sorted_pop[i] = pop[indices[i]];
    pop = sorted_pop;

}

/* Output function */

void output(vector<vector<int>> const &pop, bool file_output){

    best_ind = pop[0];
    best_fitness = fitness_func(best_ind);

    if (file_output){

        ofs_history << g << ", " << best_fitness << endl;

        if ( (g + 1) % PATH_OUT_FREQ == 0){
            for (int i=0; i<N_PTS; i++){
                best_path[i][0] = point_table[best_ind[i]][0];
                best_path[i][1] = point_table[best_ind[i]][1];
                ofs_bestpath << best_path[i][0] << ", " << best_path[i][1] << endl;
            }
            ofs_bestpath << best_path[0][0] << ", " << best_path[0][1] << endl;
        }

        if (g % DIST_OUT_FREQ == 0){
            for(int i=0; i<N_POP; i++){
                ofs_distribution << g << ", " << fitness_func(pop[i]) << endl;
            }
        }

    }

}

/* Selection */

vector<vector<int>> tournament_selection(vector<vector<int>> const &pop){

    vector<vector<int>> offspring(N_POP, vector<int>(N_PTS));
    int rand_index;
    int min_index;

    for(int i=0; i<N_POP; i++){
        min_index = N_POP;
        for (int j=0; j<N_TOURNAMENT; j++){
            rand_index = mt_engine()%N_POP;
            if (rand_index < min_index) min_index = rand_index;
        }
        offspring[i] = pop[min_index];
    }

    return offspring;
}

vector<vector<int>> roulette_selection(vector<vector<int>> const &pop, vector<double> score_vec){

    vector<vector<int>> offspring(N_POP, vector<int>(N_PTS));

    double sum = std::accumulate(score_vec.begin(), score_vec.end(), 0.0);

    for (int i=0; i<N_POP; i++){
        double rnd = get_rand_range_dbl(0.0, sum);
        for (int j=0; j<N_POP; j++){
            rnd -= score_vec[j];
            if (rnd < 0.0){
                offspring[i] = pop[j];
                break;
            }
        }
    }

    return offspring;
}

vector<vector<int>> elitist_selection(vector<vector<int>> const &pop){
    vector<vector<int>> elites(N_ELITES, vector<int>(N_PTS));
    for (int i=0; i<N_ELITES; i++) elites[i] = pop[i];
    return elites;
}

/* Representation converter */

void path_to_order(vector<int> &path){
    vector<int> indices(N_PTS);

    std::iota(indices.begin(), indices.end(), 0);
    vector<int> order = {};

    for (int i=0; i<N_PTS; i++) {
        int val = findIndex(indices, path[i]);
        order.insert(order.end(), val);
        indices.erase(indices.begin() + val);
    }

    path = order;
}

void order_to_path(vector<int> &order){
    vector<int> indices(N_PTS);
    std::iota(indices.begin(), indices.end(), 0);

    vector<int> path = {};

    for (int i=0; i<N_PTS; i++) {
        int val = indices[order[i]];
        path.insert(path.end(), val);
        indices.erase(indices.begin() + order[i]);
    }

    order = path;
}


/* Crossover operators */

void oneptcx(vector<int> &parent1, vector<int> &parent2){
    vector<int> child1={};
    vector<int> child2={};

    path_to_order(parent1);
    path_to_order(parent2);

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-1);

    // Create child.
    child1.insert(child1.begin(), parent1.begin(), parent1.begin()+pt1);
    child1.insert(child1.begin()+pt1, parent2.begin()+pt1, parent2.end());
    child2.insert(child2.begin(), parent2.begin(), parent2.begin()+pt1);
    child2.insert(child2.begin()+pt1, parent1.begin()+pt1, parent1.end());

    parent1 = child1;
    parent2 = child2;

    order_to_path(parent1);
    order_to_path(parent2);

}

void ox1(vector<int> &parent1, vector<int> &parent2){
    vector<int> child1={};
    vector<int> child2={};

    // Randomly pick the two points.
    int pt1 = get_rand_range_int(0, N_PTS-2);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS-1);

    // Extract the mid part of the vector -> "portion".
    child1.insert(child1.begin(), parent1.begin()+pt1, parent1.begin()+pt2);
    child2.insert(child2.begin(), parent2.begin()+pt1, parent2.begin()+pt2);

    // Scramble and erase.
    parent1.insert(parent1.begin(), parent1.begin() + pt2, parent1.end());
    parent1.erase(parent1.begin() + N_PTS, parent1.end());
    for (int val: child2){
        parent1.erase(parent1.begin() + findIndex(parent1, val));
    }

    parent2.insert(parent2.begin(), parent2.begin() + pt2, parent2.end());
    parent2.erase(parent2.begin() + N_PTS, parent2.end());
    for (int val: child1){
        parent2.erase(parent2.begin() + findIndex(parent2, val));
    }

    // Insert the portion at the "pt3" of the vector.
    child2.insert(child2.end(), parent1.begin(), parent1.begin()+N_PTS-pt2);
    child2.insert(child2.begin(), parent1.begin()+N_PTS-pt2, parent1.end());
    child1.insert(child1.end(), parent2.begin(), parent2.begin()+N_PTS-pt2);
    child1.insert(child1.begin(), parent2.begin()+N_PTS-pt2, parent2.end());

    parent1 = child1;
    parent2 = child2;
}

/* Mutation Operators */

void ex(vector<int> &vec){

    // Random exchange
    int r1 = mt_engine() % N_PTS;
    int r2 = mt_engine() % N_PTS;
    int tmp = vec[r1];
    vec[r1] = vec[r2];
    vec[r2] = tmp;

}

void sim(vector<int> &vec){
    // Simple Inversion Mutation

    vector<int> portion={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-2);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS-1);
    int pt3 = pt1;

    // Extract the mid part of the vector -> "portion".
    portion.insert(portion.begin(), vec.begin()+pt1, vec.begin()+pt2);
    vec.erase(vec.begin()+pt1, vec.begin()+pt2);

    // Insert the portion at the "pt3" of the vector.
    std::reverse(portion.begin(), portion.end());
    vec.insert(vec.begin()+pt3, portion.begin(), portion.end());

}


void crossover(
    std::function< void(vector<int>&, vector<int>&) > cx_operator,
    vector<vector<int>> &pop){

    std::shuffle(pop.begin(), pop.end(), mt_engine);
    for(int i=0; i<N_POP/2; i++){
        cx_operator(pop[i], pop[i+N_POP/2]);
    }
}

void mutation(
    std::function<void(vector<int>&)> mut_operator,
    vector<vector<int>> &pop){

    for (int i=0; i<N_POP; i++){
        if (get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mut_operator(pop[i]);
        }
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

    // Initial population
    vector<vector<int>> pop = create_initial_pop();

    // Evaluation
    vector<double> scores(N_POP);
    for (int i=0; i<N_POP; i++) scores[i] = fitness_func(pop[i]);

    // Sort
    sort_pop(pop, scores);

    // Output
    output(pop, true);

    for (g=1; g<NGEN; g++){

        // --- Selection
        vector<vector<int>> elites = elitist_selection(pop);
        vector<vector<int>> offspring;

        if (selection_type == "tournament"){
            offspring = tournament_selection(pop);
        }else if (selection_type == "roulette"){
            offspring = roulette_selection(pop, scores);
        }else{
            cout << "Error: Selection_type invalid" << endl;
        }

        // --- Crossover & Mutation

        if (variation_type == "simple"){
            crossover(oneptcx, offspring);
            mutation(ex, offspring);
        }else if(variation_type == "improved"){
            crossover(ox1, offspring);
            mutation(sim, offspring);
        }

        // --- offspring = offspring + elites
        offspring.insert(offspring.end(), elites.begin(), elites.end());

         // --- Evaluation
        vector<double> scores(N_POP+N_ELITES);
        for (int i=0; i < N_POP + N_ELITES ; i++) scores[i] = fitness_func(offspring[i]);

        // Sort
        sort_pop(offspring, scores);

        // Update population
        pop = offspring;

        // --- Output
        output(pop, true);

        if (problem_type == "MAX"){
            cout << "Generation: " << g << " / the longest path: " << best_fitness << endl;
        }else if(problem_type == "MIN"){
            cout << "Generation: " << g << " / the shortest path: " << 1.0/best_fitness << endl;
        }

    }

    return 0;
}
