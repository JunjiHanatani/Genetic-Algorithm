#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <functional>
#include "node.h"
#include "utility.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

vector<double> (*operator_list[])(vector<double>&, vector<double>&) = {add, sub, mul, div, sin, cos};
std::string name_list[] = {"+", "-", "*", "/", "sin", "cos"};
int const NUM_OPERATORS = 6;
double const MIN_REAL=-10.0;
double const MAX_REAL=10.0;
int const MAX_TREE_DEPTH = 5;
;
double const SNIPPING_THRESH = 1e-6;
Node const EMPTY_NODE = {"", nullptr, 0.0, -1};

vector<Node> createTree(){

    int rand_int;
    int heap_size = std::pow(2, MAX_TREE_DEPTH);
    vector<Node> nodes(heap_size);
    int parent_type = 2;
    int side = 0;
    int i = 1;

    for (int depth=0; depth<MAX_TREE_DEPTH; depth++){

        int max_width = std::pow(2, depth);

        for (int width=0; width<max_width; width++){
            Node node;

            if (depth == 0){
                parent_type = 2;
                side = 0;
            }else{
                parent_type = nodes[(int) (i/2)].type;
                side = (int) (i%2);
            }

            if (parent_type==2 || (parent_type==1 && side==0) ){

                if (depth == MAX_TREE_DEPTH - 1){
                    rand_int = get_rand_range_int(0, 1);
                }else if(depth == 0){
                    rand_int = get_rand_range_dbl(2, NUM_OPERATORS+1);
                }else{
                    rand_int = get_rand_range_dbl(0, NUM_OPERATORS+1);
                }

                if(rand_int == 0){
                    node = {"x", nullptr, 0.0, 0};
                }else if(rand_int == 1){
                    double val = get_rand_range_int(MIN_REAL, MAX_REAL);
                    node = {std::to_string(val), nullptr, val, 0};
                }else if(rand_int <= 5){
                    node = {name_list[rand_int-2], operator_list[rand_int-2], 0.0, 2};
                }else if(rand_int <= NUM_OPERATORS + 1){
                    node = {name_list[rand_int-2], operator_list[rand_int-2], 0.0, 1};
                }

            }else{
                node = EMPTY_NODE;
            }
            nodes[i] = node;
            i++;

        }

    }

    pruning(nodes);
    return nodes;
}

bool operator<(const Tree& left, const Tree& right){
  return left.fitness< right.fitness ;
}

void sort_pop(vector<Tree> &pop){
    std::sort(pop.begin(), pop.end());
}

int sample_node(vector<Node> const &tree){

    int tree_size = tree.size();
    vector<int> node_list;

    for (int i=1; i<tree_size; i++){
        if (tree[i].type != -1){
            node_list.push_back(i);
        }
    }

    int rand_int = get_rand_range_int(0, node_list.size()-1);

    return node_list[rand_int];
}

void insert_subtree(vector<Node> &tree, int const insert_pt, vector<Node> const &subtree){

    int insert_pt_depth = log2(insert_pt)+1;
    int tree_depth = log2(tree.size());
    int subtree_depth = log2(subtree.size());
    int add_depth = (insert_pt_depth - 1 + subtree_depth) - tree_depth;

    // Expand the tree depth if the subtree depth exceeds the capacity of the original tree.
    while (add_depth > 0){
        tree_depth += 1;
        add_depth -= 1;
        int width = std::pow(2, tree_depth-1);
        for (int i=0; i<width; i++) tree.push_back(EMPTY_NODE);
    }

    // Insert subtree into tree.
    int i = 1;
    for (int depth=0; depth<subtree_depth; depth++){
        int max_width = std::pow(2, depth);
        for (int width=0; width<max_width; width++){
            int tree_index = insert_pt * std::pow(2, depth) + width;
            tree[tree_index] = subtree[i];
            i++;
        }
    }
}

void pruning(vector<Node> &tree){

    int tree_depth = std::log2(tree.size());
    bool row_is_empty = true;

    for (int depth=tree_depth; depth>0; depth--){

        // Check if all entries of the row are empty.
        int max_width = std::pow(2, depth-1);
        for (int width=0; width<max_width; width++){
            int index = std::pow(2, depth-1) + width;

            // If not empty, break.
            if (tree[index].type!=-1){
                row_is_empty = false;
                break;
            }

        }

        if (row_is_empty){
            // If empty, erase all entries in the row.
            tree.erase(tree.end() - max_width, tree.end());
        }else{
            // If not empty, break.
            break;
        }
    }

}

vector<Node> pop_subtree(vector<Node> &tree, int pop_pt, bool pop){

    int pop_pt_depth = std::log2(pop_pt)+1;
    int tree_depth = std::log2(tree.size());
    int subtree_depth = tree_depth - pop_pt_depth + 1;

    vector<Node> subtree={EMPTY_NODE};
    for (int depth=1; depth<=subtree_depth; depth++){
        int max_width = std::pow(2, depth-1);
        for (int width=0; width<max_width; width++){
            int index = pop_pt * std::pow(2, depth-1) + width;
            subtree.push_back(tree[index]);

            if (pop){tree[index] = EMPTY_NODE;}
        }
    }

    if (pop) {pruning(tree);}

    return subtree;
}

vector<double> add(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x+y;});
    return v3;
};

vector<double> sub(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x-y;});
    return v3;
};

vector<double> mul(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x*y;});
    return v3;
};

vector<double> div(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x/y;});
    return v3;
};

vector<double> sin(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v3.begin(),
                   [](double x){return sin(x);});
    return v3;
};

vector<double> cos(vector<double> &v1, vector<double> &v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v3.begin(),
                   [](double x){return cos(x);});
    return v3;
};
