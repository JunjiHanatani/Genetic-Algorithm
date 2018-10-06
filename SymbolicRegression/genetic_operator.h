
#ifndef _CALC_H_
#define _CALC_H_

#include <vector>
#include "modelSymbolicRegression.h"

using std::vector;



vector<Tree> tournament_selection(vector<Tree>);

vector<Tree> roulette_selection(vector<Tree>);

vector<Tree> elitist_selection(vector<Tree>);

void crossover(std::function< void(Tree&, Tree&) > , vector<Tree>);

void mutation(std::function< void(Tree&) > , vector<Tree>);

void sort_pop(vector<Tree>);

#endif

