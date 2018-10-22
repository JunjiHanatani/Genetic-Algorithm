#ifndef _EVALUATION_H_
#define _EVALUATION_H_

#include <vector>
#include <string>
#include "node.h"

using std::vector;
using std::string;

string showFunction(Tree, bool);
vector<double> calcFunction(vector<Node> &, vector<double>);
void evaluate(vector<Tree> &);
double calcDiversity(vector<Tree>);
double calcDistance(vector<Node>const&, vector<Node>const&);
void fitnessSharing(vector<Tree> &);
#endif
