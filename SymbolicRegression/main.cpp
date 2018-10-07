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

        //auto time0 = std::chrono::system_clock::now();
        //cout<<"Selection"<<endl;
        // --- Selection
        elites = elitist_selection(pop);
        offspring = tournament_selection(pop);

        // --- Crossover
        //auto time1 = std::chrono::system_clock::now();
        crossover(oneptcx, offspring);


        // --- Mutation
        //cout<<"Mutation "<<endl;
        //auto time2 = std::chrono::system_clock::now();
        mutation(pm, offspring);

        // --- offspring + elites
        offspring.insert(offspring.end(), elites.begin(), elites.end());

        // --- Evaluation
        //cout<<"Evaluation"<<endl;
        //auto time3 = std::chrono::system_clock::now();
        evaluate(offspring);
        sort_pop(offspring);

        // Update population
        pop = offspring;

        // --- Output
        output(pop, true);
        //auto time4 = std::chrono::system_clock::now();
        /*
        std::chrono::duration<double> elapsed_seconds0 = time1-time0;
        std::chrono::duration<double> elapsed_seconds1 = time2-time1;
        std::chrono::duration<double> elapsed_seconds2 = time3-time2;
        std::chrono::duration<double> elapsed_seconds3 = time4-time3;

        cout << "Selection : " << elapsed_seconds0.count() << "s\n";
        cout << "Crossover : " << elapsed_seconds1.count() << "s\n";
        cout << "Mutation  : " << elapsed_seconds2.count() << "s\n";
        cout << "Evaluation: " << elapsed_seconds3.count() << "s\n";
        */
    }

    return 0;

}
