#ifndef _CUBE_MODEL_H_
#define _CUBE_MODEL_H_

#include <vector>
#include <chrono>
using std::vector;

// Masses
struct Mass{
  double m;
  vector<double> p;
  vector<double> v;
  vector<double> a;
};

// Springs
struct Spring{
  double k;
  double l0;
  vector<int> masses;
};

void InitializeCube(void);
void PhysicsEngine(void);

extern const int N_MASS;
extern const int N_SPRING;
extern const double nodeRadius;
extern const double ga;
extern Mass mass[];
extern Spring spring[];
extern double totalPEk;
extern double totalPEc;
extern double totalPEg;
extern double totalKE;
extern double totalEE;
extern double totalEnergy;
extern std::chrono::duration<double> elapsed_seconds;

#endif
