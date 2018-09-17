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

using std::cout;
using std::endl;
using std::vector;

// MAP
const int N_PTS = 500;
int point_table[N_PTS][2];

// GA parameter
const int NGEN = 100000;
const int N_POP = 100;
const int N_ELITES = 10;

const double P_EX = 0.0;
const double P_DM = 0.0;
const double P_IVM = 0.0;
const double P_SIM = 0.3;
const double P_ISM = 0.0;
const double P_SM = 0.0;

int g;

const int N_POS = 5;
const int N_TOURNAMENT = 3;

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

double fitness_func(vector<int> path){

    double distance = 0.0;
    for(int i=0; i<N_PTS; i++){

        int index0 = path[i];
        int index1 = path[i + 1];
        if(i==N_PTS-1){
            index1 = path[0];
        }

        int x0 = point_table[index0][0];
        int y0 = point_table[index0][1];
        int x1 = point_table[index1][0];
        int y1 = point_table[index1][1];
        distance = distance +
                   std::sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
    }
    return distance;
}

vector<vector<int>> sort_pop(vector<vector<int>> pop, vector<double> scores){
    int num = pop.size();
    vector<vector<int>> sorted_pop(num, vector<int>(N_PTS));
    vector<int> indices(num);

    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
        [&](int i, int j){return scores[i] < scores[j];}
    );

    for (int i=0; i<num; i++){
        sorted_pop[i] = pop[indices[i]];
    }

    return sorted_pop;
}

void output(vector<vector<int>> pop, bool file_output){
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

bool comp(vector<int>& lh, vector<int>& rh) {
    return fitness_func(lh) < fitness_func(rh);
}

vector<vector<int>> elite_reservation(vector<vector<int>> pop){
    vector<vector<int>> elites(N_ELITES, vector<int>(N_PTS));
    for (int i=0; i<N_ELITES; i++) elites[i] = pop[i];
    return elites;
}

vector<vector<int>> tournament_selection(vector<vector<int>> pop, int tournament_size){

    vector<vector<int>> offspring(N_POP, vector<int>(N_PTS));
    int rand_index;
    int min_index;

    for(int i=0; i<N_POP; i++){
        min_index = N_POP;
        for (int j=0; j<tournament_size; j++){
            rand_index = mt_engine()%N_POP;
            if (rand_index < min_index) min_index = rand_index;
        }
        offspring[i] = pop[min_index];
    }

    return offspring;
}

vector<vector<int>> roulette_selection(vector<vector<int>> pop, vector<double> score_vec){

    vector<vector<int>> offspring(N_POP, vector<int>(N_PTS));

    for (int i=0; i<N_POP; i++) {
        score_vec[i] = 1.0/score_vec[i];
    }

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

vector<vector<int>> pmx(vector<int> parent1, vector<int> parent2){
    vector<vector<int>> children(2, vector<int>(N_PTS));
    vector<int> child1={};
    vector<int> child2={};
    vector<int> portion1={};
    vector<int> portion2={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-2);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS-1);

    // Extract the mid part of the vector -> "portion".
    portion1.insert(portion1.begin(), parent1.begin()+pt1, parent1.begin()+pt2);
    portion2.insert(portion2.begin(), parent2.begin()+pt1, parent2.begin()+pt2);

    // Create child.
    child1.insert(child1.begin(), parent1.begin(), parent1.begin()+pt1);
    child1.insert(child1.begin()+pt1, portion2.begin(), portion2.end());
    child1.insert(child1.begin()+pt2, parent1.begin()+pt2, parent1.end());
    child2.insert(child2.begin(), parent2.begin(), parent2.begin()+pt1);
    child2.insert(child2.begin()+pt1, portion1.begin(), portion1.end());
    child2.insert(child2.begin()+pt2, parent2.begin()+pt2, parent2.end());

    // Partially-mapped operation
    for(int i=0; i<N_PTS; i++){
        if(i <pt1 || pt2 <= i){
            int gene = child1[i];
            int index = findIndex(portion2, gene);
            while (index != -1){
                gene = portion1[index];
                index = findIndex(portion2, gene);
            }
            child1[i] = gene;

            gene = child2[i];
            index = findIndex(portion1, gene);
            while (index != -1){
                gene = portion2[index];
                index = findIndex(portion1, gene);
            }
            child2[i] = gene;
        }
    }

    children[0] = child1;
    children[1] = child2;

    return children;
}

vector<vector<int>> pos(vector<int> parent1, vector<int> parent2){
    vector<vector<int>> children(2, vector<int>(N_PTS));
    vector<int> child1(N_PTS);
    vector<int> child2(N_PTS);

    vector<int> part1;
    vector<int> part2;
    vector<int> pts={};

    int pt;

    // Randomly select N points.
    do {
        pt = get_rand_range_int(0, N_PTS - 1);
        if (findIndex(pts, pt) == -1){
            pts.push_back(pt);
            //size = pts.size();
        }
    }while (pts.size() < 3);

    // Child1
    part1 = {}; part2 = parent1;
    for (int x: pts){
        part1.push_back(parent2[x]);
        int pop_index = findIndex(part2, parent2[x]);
        part2.erase(part2.begin() + pop_index);
    }
    int i1=0; int i2=0;
    for (int i=0; i<N_PTS; i++){
        if (findIndex(pts, i)!=-1){
            child1[i] = part1[i1];
            i1++;
        }else{
            child1[i] = part2[i2];
            i2++;
        }
    }

    // Child2
    part1 = {}; part2 = parent2;
    for (int x: pts){
        part1.push_back(parent1[x]);
        int pop_index = findIndex(part2, parent1[x]);
        part2.erase(part2.begin() + pop_index);
    }

    i1=0; i2=0;
    for (int i=0; i<N_PTS; i++){
        if (findIndex(pts, i)!=-1){
            child2[i] = part1[i1];
            i1++;
        }else{
            child2[i] = part2[i2];
            i2++;
        }
    }

    children[0] = child1;
    children[1] = child2;

    return children;
}

vector<vector<int>> ox1(vector<int> parent1, vector<int> parent2){
    vector<vector<int>> children(2, vector<int>(N_PTS));
    vector<int> child1={};
    vector<int> child2={};

    // Randomly pick the middle part of the vector.
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

    children[0] = child1;
    children[1] = child2;

    return children;
}

vector<vector<int>> crossover(vector<vector<int>> pop){
    vector<vector<int>> offspring(N_POP, vector<int>(N_PTS));
    vector<vector<int>> children(2, vector<int>(N_PTS));
    int rand_index1, rand_index2;

    for(int i=0; i<N_POP/2; i++){

        rand_index1 = mt_engine()%N_POP;
        rand_index2 = mt_engine()%N_POP;

        // children = pmx(pop[rand_index1], pop[rand_index2]);
        // children = pos(pop[rand_index1], pop[rand_index2]);
        children = ox1(pop[rand_index1], pop[rand_index2]);

        offspring[i] = children[0];
        offspring[i+N_POP/2] = children[1];
    }

    return offspring;
}

vector<int> ex(vector<int> vec){

    // Random exchange
    int r1 = mt_engine() % N_PTS;
    int r2 = mt_engine() % N_PTS;
    int tmp = vec[r1];
    vec[r1] = vec[r2];
    vec[r2] = tmp;

    return vec;
}

vector<int> dm(vector<int> vec){
    // Displacement Mutation
    vector<int> portion={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-1);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS);
    int pt3 = get_rand_range_int(0, N_PTS - (pt2 - pt1));

    // Extract the mid part of the vector -> "portion".
    portion.insert(portion.begin(), vec.begin()+pt1, vec.begin()+pt2);
    vec.erase(vec.begin()+pt1, vec.begin()+pt2);

    // Insert the portion at the "pt3" of the vector.
    vec.insert(vec.begin()+pt3, portion.begin(), portion.end());


    return vec;
}

vector<int> ivm(vector<int> vec){
    // Inversion Mutation

    vector<int> portion={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-2);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS-1);
    int pt3 = get_rand_range_int(0, N_PTS - (pt2 - pt1));

    // Extract the mid part of the vector -> "portion".
    portion.insert(portion.begin(), vec.begin()+pt1, vec.begin()+pt2);
    vec.erase(vec.begin()+pt1, vec.begin()+pt2);

    // Insert the portion at the "pt3" of the vector.
    std::reverse(portion.begin(), portion.end());
    vec.insert(vec.begin()+pt3, portion.begin(), portion.end());

    return vec;
}

vector<int> sim(vector<int> vec){
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

    return vec;
}

vector<int> ism(vector<int> vec){
    // Insertion Mutation

    vector<int> portion={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS - 1);
    int pt2 = pt1 + 1;
    int pt3 = get_rand_range_int(0, N_PTS - 1);

    // Extract the mid part of the vector -> "portion".
    portion.insert(portion.begin(), vec.begin()+pt1, vec.begin()+pt2);
    vec.erase(vec.begin()+pt1);

    // Insert the portion at the "pt3" of the vector.
    vec.insert(vec.begin()+pt3, portion.begin(), portion.end());

    return vec;
}

vector<int> sm(vector<int> vec){
    // Scramble Mutation

    vector<int> portion={};

    // Randomly pick the middle part of the vector.
    int pt1 = get_rand_range_int(0, N_PTS-2);
    int pt2 = get_rand_range_int(pt1 + 1, N_PTS-1);
    int pt3 = pt1;
    int pt4 = get_rand_range_int(1, pt2-pt1);

    // Extract the mid part of the vector -> "portion".
    portion.insert(portion.begin(), vec.begin()+pt1, vec.begin()+pt2);
    vec.erase(vec.begin()+pt1, vec.begin()+pt2);

    // Scramble.
    portion.insert(portion.begin(), portion.begin() + pt4, portion.end());
    portion.erase(portion.begin()+(pt2-pt1), portion.end());

    // Insert the portion at the "pt3" of the vector.
    vec.insert(vec.begin()+pt3, portion.begin(), portion.end());

    return vec;
}

vector<vector<int>> mutation(vector<vector<int>> pop){

    for (int i=0; i<N_POP; i++){

        for (int m=0; m<100; m++){
            double prob = P_EX * std::pow((1.0 - P_EX), m);
            if (get_rand_range_dbl(0.0, 1.0) > prob) break;
            pop[i] = ex(pop[i]);
        }

        if (get_rand_range_dbl(0.0, 1.0) < P_DM){
            pop[i] = dm(pop[i]);

        }else if (get_rand_range_dbl(0.0, 1.0) < P_ISM){
            pop[i] = ism(pop[i]);

        }else if (get_rand_range_dbl(0.0, 1.0) < P_IVM){
            pop[i] = ivm(pop[i]);

        }else if (get_rand_range_dbl(0.0, 1.0) < P_SIM){
            pop[i] = sim(pop[i]);

        }else if (get_rand_range_dbl(0.0, 1.0) < P_SM){
            pop[i] = sm(pop[i]);
        }

        }

    return pop;
}

int main()
{

    // Set problem

    for(int i=0; i<N_PTS; i++){
        for(int j=0; j<2; j++){
            //random value
            int r = rand()%100;
            point_table[i][j] = r;
        }
    }

    // Initial
    g = 0;
    vector<int> glist={g};
    vector<vector<int>> pop = create_initial_pop();

    // Evaluation
    vector<double> scores(N_POP);
    for (int i=0; i<N_POP; i++) scores[i] = fitness_func(pop[i]);

    // Sorting
    pop = sort_pop(pop, scores);
    //std::sort(pop.begin(), pop.end(), comp);

    // Output
    output(pop, true);

    for (g=1; g<NGEN; g++){

        // Selection
        clock_t t0 = clock();
        vector<vector<int>> elites = elite_reservation(pop);
        // vector<vector<int>> offspring = tournament_selection(pop, N_TOURNAMENT);
        vector<vector<int>> offspring = roulette_selection(pop, scores);

        // Crossover
        clock_t t1 = clock();
        offspring = crossover(offspring);

        // Mutation
        clock_t t2 = clock();
        offspring = mutation(offspring);

        // Offspring = Offspring + Elites
        clock_t t3 = clock();
        offspring.insert(offspring.end(), elites.begin(), elites.end());

        // Evaluation
        vector<double> scores(N_POP+N_ELITES);
        for (int i=0; i < N_POP + N_ELITES ; i++) scores[i] = fitness_func(offspring[i]);

        // Sorting
        offspring = sort_pop(offspring, scores);

        // Update
        pop = offspring;

        // Output
        output(pop, true);
        glist.push_back(g);
        clock_t t4 = clock();
        cout << "Generation: " << g << " / Best Fitness: " << best_fitness << endl;
        cout << "Selection = " << (double)(t1 - t0) / CLOCKS_PER_SEC << "sec.\n";
        cout << "Crossover = " << (double)(t2 - t1) / CLOCKS_PER_SEC << "sec.\n";
        cout << "Mutation  = " << (double)(t3 - t2) / CLOCKS_PER_SEC << "sec.\n";
        cout << "PostProc. = " << (double)(t4 - t3) / CLOCKS_PER_SEC << "sec.\n";
    }
    /*
    #include "matplotlibcpp.h"
    namespace plt = matplotlibcpp;

    vector<double> pt_x = {};
    vector<double> pt_y = {};
    for (int i=0; i<N_PTS; i++) pt_x.push_back(point_table[i][0]);
    for (int i=0; i<N_PTS; i++) pt_y.push_back(point_table[i][1]);

    vector<double> best_path_x = {};
    vector<double> best_path_y = {};
    for (int i=0; i<N_PTS; i++) best_path_x.push_back(best_path[i][0]);
    for (int i=0; i<N_PTS; i++) best_path_y.push_back(best_path[i][1]);

    plt::plot(glist, history);
    plt::show();

    plt::plot(pt_x, pt_y);
    plt::plot(best_path_x, best_path_y);
    plt::show();
    */

    /*
    /////////////////////////////////////////////////////////////
    // Output
    cout << "Point table: " << endl;
    for(int i=0; i<N_PTS; i++){
            cout << "  "
                 << point_table[i][0] <<", "
                 << point_table[i][1] <<endl;
    }

    getchar();

    cout << "Generation: " << g << endl;
    for(int i=0; i<N_POP; i++){
        cout << "  [";
        for (int& x: pop[i]) cout << " " << x;
        cout << "] ";
        cout << "SCORE=" << fitness_func(pop[i]) << endl;
    }

    getchar();
    cout << " --- Evaluation ---" << endl;
    cout << "Best Fitness:" << endl;
    cout << best_fitness << endl;
    cout << '\n';

    cout << "Best Path(Index): " << endl;
    for (int& x:best_ind) cout<< " " << x;
    cout << endl;
    cout << '\n';

    cout << "Best Path(Points): " << endl;
    for(int i=0; i<N_PTS; i++){
            cout << "  "
                 << point_table[i][0] <<", "
                 << point_table[i][1] <<endl;
    }
    cout << '\n';

    getchar();

    cout << " --- Selection ---" << endl;
    cout << "Elites: " << endl;
    for(int i=0; i<N_ELITES; i++){
        cout << "  [";
        for (int& x: elites[i]) cout << " " << x;
        cout << "] ";
        cout << "SCORE=" << fitness_func(elites[i]) << endl;
    }

    cout << "Selected Offspring: " << endl;
    for(int i=0; i<N_POP; i++){
        cout << "  [";
        for (int& x: offspring[i]) cout << " " << x;
        cout << "] ";
        cout << "SCORE=" << fitness_func(offspring[i]) << endl;
    }

    getchar();

    cout << " --- For analysis ---" << endl;
    cout <<"Fitness History: " << endl;
    for(vector<double> v: history) cout << v[0] << ": " << v[1] << endl;
    cout << '\n';

    cout <<"Distribution: " << endl;
    for(vector<double> v: distribution) cout << v[0] << ": " << v[1] << endl;
    cout << '\n';
    */
    /*
    vector<int> vec1 = {1, 2, 3, 4, 5, 6, 7, 8};
    vector<int> vec2 = {2, 4, 6, 8, 7, 5, 3, 1};
    vector<vector<int>> children = ox1(vec1, vec2);
    for (int x: children[0]) cout << " " << x;
    cout <<endl;
    for (int x: children[1]) cout << " " << x;
    cout <<endl;
    */
}
