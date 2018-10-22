#ifndef _NODE_H_
#define _NODE_H_

using std::vector;
using std::string;

struct Node{
    string name;
    vector<double> (*func)(vector<double>&, vector<double>&);
    double value;
    int type;
};

struct Tree{
    vector<Node> nodes;
    double fitness;
    double mae;
};

vector<double> add(vector<double> &, vector<double> &);
vector<double> sub(vector<double> &, vector<double> &);
vector<double> mul(vector<double> &, vector<double> &);
vector<double> div(vector<double> &, vector<double> &);
vector<double> sin(vector<double> &, vector<double> &);
vector<double> cos(vector<double> &, vector<double> &);

void sort_pop(vector<Tree>&);
int sample_node(vector<Node> const &);
void insert_subtree(vector<Node> &, int const, vector<Node> const &);
vector<Node> pop_subtree(vector<Node> &, int, bool);
void pruning(vector<Node> &);
vector<Node> createTree();


extern vector<double> (*operator_list[])(vector<double>&, vector<double>&);
extern string name_list[];
extern int const NUM_OPERATORS;
extern double const MIN_REAL;
extern double const MAX_REAL;
extern Node const EMPTY_NODE;
extern int const MAX_TREE_DEPTH;
extern double const SNIPPING_THRESH;

#endif
