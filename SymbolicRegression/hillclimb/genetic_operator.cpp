#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "genetic_operator.h"
#include "node.h"
#include "utility.h"

int const NGEN = 1000;
int const N_POP = 1;
int const N_TOURNAMENT = 5;
int const N_ELITES = 5;
double const P_MUT = 1.0;

////////////// Initial population ///////////////

vector<Tree> create_initial_pop(
    std::function<vector<Node>()> createTree){

    vector<Tree> pop(N_POP);

    //Create population
    for (int i=0; i<N_POP; i++){
        Tree tree;
        tree.nodes = createTree();
        tree.fitness = 0.0;
        pop[i] = tree;
    }

    return pop;
}


////////////// Selection ///////////////

vector<Tree> tournament_selection(vector<Tree> const &pop){

    vector<Tree> offspring(N_POP);
    int rand_index;
    int min_index;

    for(int i=0; i<N_POP; i++){
        min_index = N_POP;
        for (int j=0; j<N_TOURNAMENT; j++){
            rand_index = get_rand_range_int(0, N_POP-1);
            if (rand_index < min_index) min_index = rand_index;
        }
        offspring[i] = pop[min_index];
    }

    return offspring;
}

vector<Tree> roulette_selection(vector<Tree> pop){

    vector<Tree> offspring;
    vector<double> rand_list(N_POP);

    // Calculate sum of the fitness over all population.
    double sum_fitness = std::accumulate(pop.begin(), pop.end(), 0.0,
                     [](double sum, Tree& tree ){ return sum+1.0/tree.fitness; } );
    // Generate random list.
    for (int i=0; i<N_POP; i++){
        rand_list[i] = get_rand_range_dbl(0.0, sum_fitness);
    }
    // Sort random_list.
    std::sort(rand_list.begin(), rand_list.end());
    std::reverse(rand_list.begin(), rand_list.end());

    double thresh = 0.0;
    for (int i=0; i<N_POP; i++){
        thresh += 1.0/pop[i].fitness;
        while(offspring.size()<N_POP){
            double rand = rand_list.back();
            if (rand<thresh){
                offspring.push_back(pop[i]);
                rand_list.pop_back();
            }else{
                break;
            }
        }
    }

    return offspring;
}

vector<Tree> elitist_selection(vector<Tree> const &pop){
    vector<Tree> elites(N_ELITES);
    for (int i=0; i<N_ELITES; i++) elites[i] = pop[i];
    return elites;
}


////////////// Mutation Operators //////////

void mutation(
    std::function<void(Tree&)> mut_operator,
    vector<Tree> &pop){

    for (int i=0; i<N_POP; i++){
        if (get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mut_operator(pop[i]);
        }
    }
}

// --- Point Mutation
void pm(Tree &tree){

    int index = sample_node(tree.nodes);
    Node target_node = tree.nodes[index];

    if (target_node.type==0){
        if (target_node.name != "x"){
            double val = get_rand_range_dbl(MIN_REAL, MAX_REAL);
            //vector<double> vals(N_PTS, val);
            target_node.value = val;
            target_node.name = std::to_string(val);
        }
    }else if(target_node.type==2){
        int randint = get_rand_range_int(0, 3);
        target_node.func = operator_list[randint];
        target_node.name = name_list[randint];

    }else if(target_node.type==1){
        int randint = get_rand_range_int(4, NUM_OPERATORS-1);
        target_node.func = operator_list[randint];
        target_node.name = name_list[randint];
    }

    tree.nodes[index] = target_node;
    tree.fitness = -1.0;
}

// --- Subtree Mutation
void sm(Tree &tree){

    int index = sample_node(tree.nodes);

    pop_subtree(tree.nodes, index, true);
    vector<Node> subtree = createTree();
    insert_subtree(tree.nodes, index, subtree);

    tree.fitness = -1.0;
}

////////////// Crossover operators /////////

void crossover(
    std::function< void(Tree&, Tree&) > cx_operator,
    vector<Tree> &pop){

    std::shuffle(pop.begin(), pop.end(), mt_engine);
    for(int i=0; i<N_POP/2; i++){
        cx_operator(pop[i], pop[i+N_POP/2]);
    }

}

// --- Single-point crossover
void oneptcx(Tree &tree1, Tree &tree2){
    vector<Node> nodes1 = tree1.nodes;
    vector<Node> nodes2 = tree2.nodes;

    // Get cut point randomly.
    int cut_pt1 = sample_node(nodes1);
    int cut_pt2 = sample_node(nodes2);

    // Get subtrees.
    vector<Node> subtree1 = pop_subtree(nodes1, cut_pt1, true);
    vector<Node> subtree2 = pop_subtree(nodes2, cut_pt2, true);

    // Swap.
    insert_subtree(nodes1, cut_pt1, subtree2);
    insert_subtree(nodes2, cut_pt2, subtree1);

    // Update.
    tree1 = {nodes1, -1.0};
    tree2 = {nodes2, -1.0};
}
