#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "evaluation.h"
#include "genetic_operator.h"
#include "io.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

string showFunction(Tree tree, bool verbose){
    vector<Node> nodes = tree.nodes;
    int imax = nodes.size()-1;

    if (imax ==0){
        return "No nodes.";
    }

    for (int i=imax; i>0; i--){

        if(nodes[i].type==1){
            nodes[i].name = nodes[i].name + "(" + nodes[2*i].name + ")";
        }else if(nodes[i].type==2){
            nodes[i].name = "(" + nodes[2*i].name + nodes[i].name + nodes[2*i+1].name + ")";
        }

    }
    if (verbose){
        cout << " f(x) = " << nodes[1].name
             << " : fitness=" << tree.fitness
             << " : Tree size= " << tree.nodes.size()
             << endl;
    }
    return nodes[1].name;
}

vector<double> calcFunction(vector<Node> &tree, vector<double>x_vec){

    int x_size = x_vec.size();
    int tree_size = tree.size();
    int last_index = tree_size-1;

    vector<vector<double>> tree_for_calc(tree_size, vector<double>(x_size));
    for (int i=0; i< tree_size; i++){
        double val = tree[i].value;
        vector<double> vals(x_size, val);
        tree_for_calc[i] = vals;
    }

    vector<double> a(x_size); vector<double> b(x_size);
    vector<double> ans(x_size);

    if (tree[1].name=="x") {return x_vec;}

    for (int i=last_index; i>0; i--){
        //if (g==37){cout << i << ":"<< tree[i].name << endl;}

        if(tree[i].type==1){

            if(tree[2*i].name=="x"){
                a = x_vec;
            }else{
                a = tree_for_calc[2*i];
            }

            ans = tree[i].func(a, a);

        }else if(tree[i].type==2){

            if(tree[2*i].name=="x"){
                a = x_vec;
            }else{
                a = tree_for_calc[2*i];
            }

            if(tree[2*i+1].name=="x"){
                b = x_vec;
            }else{
                b = tree_for_calc[2*i+1];
            }

            ans =  tree[i].func(a, b);

        }else{
            continue;
        }

        tree_for_calc[i] = ans;

        // Snipping
        double ans_max = *std::max_element(ans.begin(), ans.end());
        double ans_min = *std::min_element(ans.begin(), ans.end());

        if(((ans_max - ans_min) < SNIPPING_THRESH) || (log2(i) + 1 >= MAX_TREE_DEPTH) ){
            double ans_ave = std::accumulate(ans.begin(), ans.end(), 0.0)/x_size;
            tree[i] = {std::to_string(ans_ave), nullptr, ans_ave, 0};
            vector<double> vals(x_size, ans_ave);
            tree_for_calc[i] = vals;

            // Delete:
            vector<int> check_list = {2*i, 2*i+1};
            vector<int> next_list;
            while(check_list[0]<tree_size){
                for (int index: check_list){
                    tree[index] = EMPTY_NODE;
                    next_list.push_back(2*index);
                    next_list.push_back(2*index+1);
                }
                check_list = next_list;
                next_list.clear();
            }
        }
    }

    pruning(tree);

    return tree_for_calc[1];
}

void evaluate(vector<Tree> &pop){
    double err=0.0;
    vector<double> x_vec(N_PTS);
    vector<double> y_vec(N_PTS);

    for (int i=0; i<N_PTS; i++){
        x_vec[i] = point_table[i][0];
        y_vec[i] = point_table[i][1];
    }

    for (int i=0; i<N_POP; i++) {
        //if (g==37){
        //    cout << "------------- " << i << " ------------ " << endl;
        //}

        vector<double> y_vec_estimate = calcFunction(pop[i].nodes, x_vec);

        // Calc. root mean square
        err = 0.0;
        for (int j=0; j<N_PTS; j++){
            err += (y_vec_estimate[j] - y_vec[j]) * (y_vec_estimate[j] - y_vec[j]);
        }
        err = sqrt(err);

        // Set fitness.
        if (std::isnan(err)){
            pop[i].fitness = std::numeric_limits<double>::infinity();
        }else if(pop[i].nodes.size()<=2){
            pop[i].fitness = std::numeric_limits<double>::infinity();
        }else{
            pop[i].fitness = err;
        }
        //showFunction(pop[i], true);
    }
}

