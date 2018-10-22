#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>
#include <chrono>
#include "utility.h"
#include "node.h"
#include "evaluation.h"
#include "genetic_operator.h"
#include "io.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

int g = 0;
int MAX_EVALS = 1500000;

/* One step hill-climb */
void hillclimb(Tree &best_tree, bool &climbing_is_done){

    bool move_made = false;
    vector<Tree> pop(1);
    Tree tree;

    for (int i=0; i<100; i++){

        pop[0] = best_tree;

        // Exchange two cities.
        mutation(pm, pop);

        // Evaluate the new vector.
        evaluate(pop);
        g += 1;

        // if the best vector updates, move on next climb.
        if (pop[0].fitness < best_tree.fitness){
            best_tree = pop[0];
            move_made = true;
            break;
        }

        // if num_evals reaches the max number, the calculation ends.
        if (g >= MAX_EVALS){
            climbing_is_done = true;
            break;
        }

        if (move_made == true) break;
    }

    if (move_made == false){
        climbing_is_done = true;
    }

}


int main(){
    vector<Tree> pop(1);
    Tree best_tree;
    Tree tree;

    pop = create_initial_pop(createTree);
    evaluate(pop);
    best_tree = pop[0];

    // Read data.
    read_data();

    // Random restart hill-climbing
    bool climbing_is_done;

    while (g < MAX_EVALS){

        // Random restart
        pop = create_initial_pop(createTree);
        evaluate(pop);
        climbing_is_done = false;
        tree = pop[0];

        // Hill-climbing start
        while (climbing_is_done==false){

            // One step climbing
            hillclimb(tree, climbing_is_done);

            if(best_tree.fitness > tree.fitness){
                best_tree = tree;
            }

            pop[0] = best_tree;
            output(pop, true);
        }

    }

    return 0;
}
