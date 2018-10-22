#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "evaluation.h"
#include "genetic_operator.h"
#include "node.h"
#include "utility.h"

int const NGEN = 2000;
int const N_POP = 500;
int const N_TOURNAMENT = 5;
int const N_ELITES = 5;
double const P_MUT = 0.3;

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

    int n_offspring = N_POP - N_ELITES;
    vector<Tree> offspring(n_offspring);
    int rand_index;
    int min_index;

    for(int i=0; i<n_offspring; i++){
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
        //if (get_rand_range_dbl(0.0, 1.0) < P_MUT){
        //    mut_operator(pop[i]);
        //}

        for(int j=0; j<5; j++){
            double p = pow(1.0 - P_MUT, j) * P_MUT;
            if (get_rand_range_dbl(0.0, 1.0) < p){
                mut_operator(pop[i]);
            }
        }

    }
}

// --- Point Mutation
void pm(Tree &tree){

    int index = sample_node(tree.nodes);
    Node target_node = tree.nodes[index];

    if (target_node.type==0){
        if (get_rand_range_dbl(0, 1) < 0.3){
            target_node.value = 0.0;
            target_node.name = "x";
        }else{
            double val = get_rand_range_int(MIN_REAL, MAX_REAL);
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

// --- Single-point crossover with Crowding
void crowding_cx(Tree &tree1, Tree &tree2){
    vector<Node> nodes1 = tree1.nodes;
    vector<Node> nodes2 = tree2.nodes;
    vector<Tree> pop(4);
    vector<Node> parent1 = nodes1;
    vector<Node> parent2 = nodes2;

    // Get cut point randomly.
    int cut_pt1 = sample_node(nodes1);
    int cut_pt2 = sample_node(nodes2);

    // Get subtrees.
    vector<Node> subtree1 = pop_subtree(nodes1, cut_pt1, true);
    vector<Node> subtree2 = pop_subtree(nodes2, cut_pt2, true);

    // Swap.
    insert_subtree(nodes1, cut_pt1, subtree2);
    insert_subtree(nodes2, cut_pt2, subtree1);

    // Child
    vector<Node> child1 = nodes1;
    vector<Node> child2 = nodes2;

    pop[0] = {parent1, -1.0};
    pop[1] = {parent2, -1.0};
    pop[2] = {child1, -1.0};
    pop[3] = {child2, -1.0};

    evaluate(pop);
    if (calcDistance(parent1, child1) + calcDistance(parent2, child2) <
        calcDistance(parent1, child2) + calcDistance(parent2, child1)){

        if (pop[0].fitness < pop[2].fitness){
            tree1 = pop[0];
        }else{
            tree1 = pop[2];
        }

        if (pop[1].fitness < pop[3].fitness){
            tree2 = pop[1];
        }else{
            tree2 = pop[3];
        }

    }else{

        if (pop[0].fitness < pop[3].fitness){
            tree1 = pop[0];
        }else{
            tree1 = pop[3];
        }

        if (pop[1].fitness < pop[2].fitness){
            tree2 = pop[1];
        }else{
            tree2 = pop[2];
        }
    }

    // Update.
    //tree1 = {nodes1, -1.0};
    //tree2 = {nodes2, -1.0};

}
