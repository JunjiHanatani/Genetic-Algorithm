#ifndef GENETICOPERATORS_H
#define GENETICOPERATORS_H

#include <functional>
#include <vector>
#include "PhysicsEngine.h"
#include "CubeGenerator.h"

using std::vector;

struct Individual{
    vector<vector<double>> para;
    vector<vector<double>> trajectory;
    double fitness;
    double distance;
    double err;
    double diameter;
    int num_of_masses;
    int num_of_springs;
    int age;

    Individual():
        para(),
        trajectory(),
        fitness(),
        distance(),
        err(),
        diameter(),
        num_of_masses(),
        num_of_springs(),
        age(1)
        {}
};

extern int gen;
extern int SIZE_OF_CHROMOSOME;
extern int SIZE_OF_GENE;

/* --- Functions --- */

void alpsGA(vector<Individual>[]);
void EvolveCube(void);
void Restart(void);

// Create population
vector<Individual> createInitialPop(int);

// Evaluation
void evaluate(vector<Individual>&);
void evaluateMPI(vector<Individual>&);
void sort_pop(vector<Individual>&);
double calcFitness(double, double, double);
vector<double> calcTraveleLength(vector<vector<double>>&);
double calcDiameter(vector<Mass>&);

// Selection
vector<Individual> tournamentSelection(vector<Individual> const &, int, int);
vector<Individual> rouletteSelection(vector<Individual> &, int);
vector<Individual> elitistSelection(vector<Individual> const &, int);
vector<Individual> overageSelection(vector<Individual>&, int);
vector<Individual> agebaseElitistSelection(vector<Individual> const &, int);

// Mutation
void mutation(vector<Individual> &);
void mutNormal(Individual&);
void mutSphere(Individual&);
void mutExist(Individual&);
void mutRadius(Individual&);

// Crossover
void crossover(vector<Individual> &, vector<Individual> &);
vector<Individual> oneptcx(Individual &, const Individual &);
vector<Individual> oneptswap(Individual &, const Individual &);
vector<Individual> area_cx(Individual &, const Individual &);

void RecordLog(vector<Individual>[]);
void setRepresentation(char);
int readRestartData(vector<Individual>[], string);

#endif

