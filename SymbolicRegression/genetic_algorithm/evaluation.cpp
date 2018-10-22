#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include "evaluation.h"
#include "genetic_operator.h"
#include "io.h"
#include "utility.h"

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
             //<< " : fitness=" << tree.fitness
             //<< " : Tree size= " << tree.nodes.size()
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
            if (ans_ave > MAX_REAL){ans_ave=MAX_REAL;};
            if (ans_ave < MIN_REAL){ans_ave=MIN_REAL;};
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

    int num = pop.size();
    double mean_abs_err;
    vector<double> x_vec(N_PTS);
    vector<double> y_vec(N_PTS);
    vector<double> y_vec_estimate(N_PTS);
    vector<double> y_abs_err(N_PTS);

    for (int i=0; i<N_PTS; i++){
        x_vec[i] = point_table[i][0];
        y_vec[i] = point_table[i][1];
    }

    for (int i=0; i<num; i++) {

        vector<double> y_vec_estimate = calcFunction(pop[i].nodes, x_vec);

        // Calc. mean absolute error.
        std::transform(y_vec.begin(), y_vec.end(), y_vec_estimate.begin(), y_abs_err.begin(),
                [](double a, double b){return fabs(a-b);});
        mean_abs_err = std::accumulate(y_abs_err.begin(), y_abs_err.end(), 0.0)/(double)N_PTS;

        // Set fitness.
        if (std::isnan(mean_abs_err)){
            pop[i].fitness = std::numeric_limits<double>::infinity();
            pop[i].mae = std::numeric_limits<double>::infinity();
        }else if(pop[i].nodes.size()<=2){
            pop[i].fitness = std::numeric_limits<double>::infinity();
            pop[i].mae = std::numeric_limits<double>::infinity();
        }else{
            pop[i].fitness = mean_abs_err;
            pop[i].mae = mean_abs_err;
        }
        //showFunction(pop[i], true);
    }
    //fitnessSharing(pop);
    //for (Tree t:pop){
    //    std::cout << t.mae << " " << t.fitness << std::endl;
    //}
    //calcDiversity(pop);
}

double calcDistance(vector<Node> const &nodes1, vector<Node> const &nodes2){

    int tree_size1 = nodes1.size();
    int tree_size2 = nodes2.size();
    int tree_size_min = std::min(tree_size1, tree_size2);
    int identity_count = 0;

    for (int i=0; i< tree_size_min; i++){
        if (nodes1[i].name == nodes2[i].name){
            identity_count += 1;
        //}else if(nodes1[i].type==0 && nodes1[i].name!="x" &&
        //         nodes2[i].type==0 && nodes2[i].name!="x"){
        //    identity_count += 1;
        }
    }
    double distance = 1.0 - (double)identity_count/(double)tree_size_min;

    return distance;
}

void fitnessSharing(vector<Tree> &pop){
    double thresh = 0.5;
    double alpha = 0.5;
    int num = pop.size();
    for(int i=0; i<num; i++){
        double sum_share = 0.0;
        for(int j=0; j<num; j++){
            double distance = calcDistance(pop[i].nodes, pop[j].nodes);
            if (distance < thresh){
                sum_share += 1.0 - pow(distance/thresh, alpha);
            }
        }
        pop[i].fitness = pop[i].fitness * sum_share;
    }
}

double calcDiversity(vector<Tree> pop){

    int num = pop.size();
    double sum_distance = 0.0;
    int cnt = 0;
    for(int i=0; i<num; i++){
        for(int j=i+1; j<num; j++){
            sum_distance += calcDistance(pop[i].nodes, pop[j].nodes);
            cnt += 1;
        }
    }
    double diversity = (double)sum_distance / (double)cnt;

    return diversity;
}
