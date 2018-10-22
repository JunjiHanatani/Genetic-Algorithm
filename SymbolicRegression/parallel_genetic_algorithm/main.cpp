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

int g=0;

int main()
{
    vector<Tree> pop_total;
    vector<vector<Tree>> pop(NUM_PARALLEL);
    vector<vector<Tree>> elites(NUM_PARALLEL, vector<Tree>(N_ELITES));
    vector<vector<Tree>> offspring(NUM_PARALLEL);

    // Read data.
    read_data();

    // Initial population
    for (int i=0; i<NUM_PARALLEL; i++){
        pop[i] = create_initial_pop(createTree);
    }

    // Evaluation & sort.
    for (int i=0; i<NUM_PARALLEL; i++){
        evaluate(pop[i]);
        sort_pop(pop[i]);
    }

    // Merge & output
    pop_total = pop[0];
    for (int i=1; i<NUM_PARALLEL; i++){
        pop_total.insert(pop_total.end(), pop[i].begin(), pop[i].end());
    }
    sort_pop(pop_total);
    output(pop_total, true);


    for (g=1; g<=NGEN; g++){

        // --- Selection
        for (int i=0; i<NUM_PARALLEL; i++){
            elites[i] = elitist_selection(pop[i]);
            offspring[i] = tournament_selection(pop[i]);
            //offspring = roulette_selection(pop);
        }

        // --- Crossover
        for (int i=0; i<NUM_PARALLEL; i++){
            crossover(oneptcx, offspring[i]);
        }

        // --- Mutation
        for (int i=0; i<NUM_PARALLEL; i++){
            mutation(pm, offspring[i]);
        }

        // --- offspring + elites
        for (int i=0; i<NUM_PARALLEL; i++){
            offspring[i].insert(offspring[i].end(), elites[i].begin(), elites[i].end());
        }

        // Evaluation & sort.
        for (int i=0; i<NUM_PARALLEL; i++){
            evaluate(offspring[i]);
            sort_pop(offspring[i]);
        }

        // Update population
        for (int i=0; i<NUM_PARALLEL; i++){
            pop[i] = offspring[i];
        }

        // Merge & output
        pop_total = pop[0];
        for (int i=1; i<NUM_PARALLEL; i++){
            pop_total.insert(pop_total.end(), pop[i].begin(), pop[i].end());
        }
        sort_pop(pop_total);
        output(pop_total, true);

    }

    return 0;

}
