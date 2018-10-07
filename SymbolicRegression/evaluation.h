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

#endif
