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
    vector<Tree> pop;
    vector<Tree> elites(N_ELITES);
    vector<Tree> offspring;

    // Read data.
    read_data();

    // Initial population
    pop = create_initial_pop(createTree);

    // Evaluation & sort.
    evaluate(pop);
    sort_pop(pop);

    // Output
    output(pop, true);


    for (g=1; g<=NGEN; g++){

        // --- Selection
        elites = elitist_selection(pop);
        offspring = tournament_selection(pop);
        offspring.insert(offspring.end(), elites.begin(), elites.end());
        //offspring = roulette_selection(pop);

        // --- Crossover
        crossover(oneptcx, offspring);
        //crossover(crowding_cx, offspring);

        // --- Mutation
        mutation(pm, offspring);

        // --- offspring + elites
        offspring.insert(offspring.end(), elites.begin(), elites.end());

        // --- Evaluation
        evaluate(offspring);
        sort_pop(offspring);

        // Update population
        pop = offspring;

        // --- Output
        output(pop, true);
    }

    return 0;

}
