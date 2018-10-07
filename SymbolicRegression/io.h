#ifndef _IO_H_
#define _IO_H_

#include <vector>
#include "node.h"

using std::vector;

extern int const DOWN_SAMPLE_STEP;
extern int const N_PTS;
extern std::ifstream ifs;

// Results
extern int const DIST_OUT_FREQ;
extern int const PATH_OUT_FREQ;
extern std::ofstream ofs_history;
extern std::ofstream ofs_distribution;
extern std::ofstream ofs_bestpath;
extern vector<vector<double>> point_table;

void read_data();
void output(vector<Tree>, bool);

#endif

