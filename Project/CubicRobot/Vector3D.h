#ifndef Vector3D_H
#define Vector3D_H
#include <vector>
using std::vector;

double calcDistance(vector<double>&, vector<double>&);
double calcNorm(vector<double>&);
vector<double> calcUnitVector(vector<double>&);
vector<double> calcUnitVector2(vector<double>&, vector<double>&);
double dotProduct(vector<double>&, vector<double>&);
vector<double> calcVector(vector<double>&, double);
vector<double> add(vector<double>&, vector<double>&);
vector<double> sub(vector<double>&, vector<double>&);
vector<double> mul(vector<double>&, vector<double>&);
vector<double> div(vector<double>&, vector<double>&);
vector<double> scaling(vector<double>&, double);

#endif
