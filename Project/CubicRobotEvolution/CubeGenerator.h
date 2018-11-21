#ifndef _CUBEGENERATOR_H_
#define _CUBEGENERATOR_H_

#include <vector>
#include <string>
#include <array>
#include <GL/glut.h>

using std::vector;
using std::string;

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

extern vector<Mass> mass;
extern vector<Spring> spring;
extern vector<std::array<int, 8>> vertices;

extern int N_MASS;
extern int N_SPRING;
extern int N_CUBE;
extern int N_SYMMETRIC_PAIR;

extern const double nodeRadius;
extern vector<double> springInitialRestlength;
extern vector<double> breathe_amp;
extern vector<double> breathe_phase;
extern vector<double> breathe_offset;
extern const double breathe_omega;


/* --- Functions --- */
void InitializeCube(void);
void GenerateCube(void);
void SetBreathe(vector<vector<double>>, string);
void FindSymmetricSpring(void);
void CubeLog(void);

#endif
