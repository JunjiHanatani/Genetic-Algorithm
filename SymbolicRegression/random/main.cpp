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


int main(){
    vector<Tree> pop(1);
    Tree best_tree;
    Tree tree;

    pop = create_initial_pop(createTree);
    evaluate(pop);
    best_tree = pop[0];

    // Read data.
    read_data();


    for (g=0; g < MAX_EVALS; g++){

        // New solution
        pop = create_initial_pop(createTree);
        evaluate(pop);
        tree = pop[0];

        // Update.
        if (tree.fitness < best_tree.fitness){
            best_tree = tree;
        }

        // --- Output
        pop[0] = best_tree;
        output(pop, true);

    }

    return 0;
}
