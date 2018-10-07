#ifndef _GENETIC_OPERATOR_H_
#define _GENETIC_OPERATOR_H_

#include <vector>
#include <functional>
#include "node.h"

using std::vector;

vector<Tree> create_initial_pop(std::function<vector<Node>()>);

vector<Tree> tournament_selection(vector<Tree> const &);
vector<Tree> roulette_selection(vector<Tree> const &);
vector<Tree> elitist_selection(vector<Tree> const &);

void crossover(std::function< void(Tree&, Tree&) > , vector<Tree>&);
void oneptcx(Tree &, Tree &);

void mutation(std::function< void(Tree&) > , vector<Tree>&);
void pm(Tree &);

extern int const NGEN;
extern int const N_POP;
extern int const N_TOURNAMENT;
extern int const N_ELITES;
extern double const P_MUT;
extern int g;

#endif
