#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <functional>
#include "utility.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using std::cout;
using std::endl;
using std::vector;
using std::string;

struct Node{
    std::string name;
    double (*func)(double, double);
    double value;
    int type;
};

struct Tree{
    vector<Node> nodes;
    double fitness;
};

int const NGEN = 10;
int const N_PTS = 1000;
int const N_POP = 100;
int const N_TOURNAMENT = 3;
int const N_ELITES = 3;
double const P_MUT=0.3;
int g;

int const MAX_TREE_DEPTH = 3;
int const HEAP_SIZE = std::pow(2, MAX_TREE_DEPTH);


std::ifstream ifs("./sample/function1.csv");
vector<vector<double>> point_table(N_PTS, vector<double>(2));

// Results
int DIST_OUT_FREQ = NGEN;
int PATH_OUT_FREQ = NGEN;
std::ofstream ofs_history("./history.csv");
std::ofstream ofs_distribution("./distribution.csv");
std::ofstream ofs_bestpath("./bestpath.csv");

double XMIN = 0.0;
double XMAX = 10.0;
double DX = 0.1;


double add(double val1, double val2){return val1 + val2;}
double sub(double val1, double val2){return val1 - val2;}
double mul(double val1, double val2){return val1 * val2;}
double div(double val1, double val2){return val1 / val2;}
double sin(double val1, double val2){return sin(val1);}
double cos(double val1, double val2){return cos(val1);}
double tan(double val1, double val2){return tan(val1);}
double log(double val1, double val2){return log(val1);}
double exp(double val1, double val2){return exp(val1);}
double sqr(double val1, double val2){return sqrt(val1);}

double (*operator_list[])(double, double) =
    {add, sub, mul, div, sin, cos};

std::string name_list[] =
    {"+", "-", "*", "/", "sin", "cos"};

int NUM_OPERATORS = 6;
double MIN_REAL=-10.0;
double MAX_REAL=10.0;


// Selection

vector<Tree> tournament_selection(vector<Tree> const &pop){

    vector<Tree> offspring(N_POP);
    int rand_index;
    int min_index;

    for(int i=0; i<N_POP; i++){
        min_index = N_POP;
        for (int j=0; j<N_TOURNAMENT; j++){
            rand_index = get_rand_range_int(0, N_POP-1);
            if (rand_index < min_index) min_index = rand_index;
        }
        offspring[i] = pop[min_index];
    }

    return offspring;
}

vector<Tree> roulette_selection(vector<Tree> const &pop){

    vector<Tree> offspring(N_POP);
    vector<double> rand_list(N_POP);

    // Calculate sum of the fitness over all population.
    //double sum = std::accumulate(pop.begin(), pop.end(), 0.0,
    //                 []( int sum, Tree& tree ){ return sum+tree.fitness; } );
    double sum=10.0;
    // Generate random list.
    for (int i=0; i<N_POP; i++){
        rand_list[i] = get_rand_range_dbl(0.0, sum);
    }
    // Sort random_list.
    std::sort(rand_list.begin(), rand_list.end());
    std::reverse(rand_list.begin(), rand_list.end());

    double thresh = 0.0;
    for (int i=0; i<N_POP; i++){
        thresh += pop[i].fitness;
        while(offspring.size()<N_POP){
            double rand = rand_list.back();
            if (rand<thresh){
                offspring.push_back(pop[i]);
                rand_list.pop_back();
            }else{
                break;
            }
        }
    }

    return offspring;
}

vector<Tree> elitist_selection(vector<Tree> const &pop){
    vector<Tree> elites(N_ELITES);
    for (int i=0; i<N_ELITES; i++) elites[i] = pop[i];
    return elites;
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
    Node empty_node = {"", nullptr, 0.0, -1};

    // Expand the tree depth if the subtree depth exceeds the capacity of the original tree.
    while (add_depth > 0){
        tree_depth += 1;
        add_depth -= 1;
        int width = std::pow(2, tree_depth-1);
        for (int i=0; i<width; i++) tree.push_back(empty_node);
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

vector<Node> pop_subtree(vector<Node> &tree, int pop_pt, bool pop){

    int pop_pt_depth = std::log2(pop_pt)+1;
    int tree_depth = std::log2(tree.size());
    int subtree_depth = tree_depth - pop_pt_depth + 1;

    Node empty_node = {"", nullptr, 0.0, -1};
    vector<Node> subtree={empty_node};
    for (int depth=1; depth<=subtree_depth; depth++){
        int max_width = std::pow(2, depth-1);
        for (int width=0; width<max_width; width++){
            int index = pop_pt * std::pow(2, depth-1) + width;
            subtree.push_back(tree[index]);
            if (pop){
                tree[index] = empty_node;
            }
        }
    }


    if (pop){
        bool row_is_empty = true;
        for (int depth=tree_depth; depth>0; depth--){
            // Check if all entries of the row are empty.
            int max_width = std::pow(2, depth-1);
            for (int width=0; width<max_width; width++){
                int index = std::pow(2, depth-1) + width;
                if (tree[index].type!=-1){
                    row_is_empty = false;
                    break;
                }
            }

            if (row_is_empty){
                tree.erase(tree.end() - max_width, tree.end());
            }else{
                break;
            }
        }
    }

    return subtree;
}

// Mutation Operators


// Point Mutation
void pm(Tree &tree){

    int index = sample_node(tree.nodes);
    Node target_node = tree.nodes[index];

    if (target_node.type==0){
        if (target_node.name != "x"){
            double val = get_rand_range_dbl(MIN_REAL, MAX_REAL);
            target_node.value = val;
            target_node.name = std::to_string(val);
        }
    }else if(target_node.type==2){
        int randint = get_rand_range_int(0, 3);
        target_node.func = operator_list[randint];
        target_node.name = name_list[randint];

    }else if(target_node.type==1){
        int randint = get_rand_range_int(4, NUM_OPERATORS-1);
        target_node.func = operator_list[randint];
        target_node.name = name_list[randint];
    }

    tree.nodes[index] = target_node;
    tree.fitness = -1.0;
}


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
        cout << " f(x) = " << nodes[1].name << " / fitness=" << tree.fitness<< endl;
    }
    return nodes[1].name;
}

// Crossover operators

void oneptcx(Tree &tree1, Tree &tree2){

    vector<Node> nodes1 = tree1.nodes;
    vector<Node> nodes2 = tree2.nodes;

    // Get cut point randomly.
    int cut_pt1 = sample_node(nodes1);
    int cut_pt2 = sample_node(nodes2);

    // Get subtrees.
    vector<Node> subtree1 = pop_subtree(nodes1, cut_pt1, true);
    vector<Node> subtree2 = pop_subtree(nodes2, cut_pt2, true);

    // Swap.
    insert_subtree(nodes1, cut_pt1, subtree2);
    insert_subtree(nodes2, cut_pt2, subtree1);

    // Update.
    tree1 = {nodes1, -1.0};
    tree2 = {nodes2, -1.0};

}


void crossover(
    std::function< void(Tree&, Tree&) > cx_operator,
    vector<Tree> &pop){

    std::shuffle(pop.begin(), pop.end(), mt_engine);
    for(int i=0; i<N_POP/2; i++){
        cx_operator(pop[i], pop[i+N_POP/2]);
    }
}

void mutation(
    std::function<void(Tree&)> mut_operator,
    vector<Tree> &pop){

    for (int i=0; i<N_POP; i++){
        if (get_rand_range_dbl(0.0, 1.0) < P_MUT){
            mut_operator(pop[i]);
        }
    }

}


vector<Tree> create_initial_pop(
    std::function<vector<Node>()> createTree){

    vector<Tree> pop(N_POP);

    //Create population
    for (int i=0; i<N_POP; i++){
        Tree tree;
        tree.nodes = createTree();
        tree.fitness = 0.0;
        pop[i] = tree;
    }

    return pop;
}


vector<Node> createTree(){

    int rand_int;
    vector<Node> nodes(HEAP_SIZE);
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
                }else{
                    rand_int = get_rand_range_int(0, NUM_OPERATORS+1);
                }

                if(rand_int == 0){
                    node = {"x", nullptr, 0.0, 0};
                }else if(rand_int == 1){
                    double val = get_rand_range_dbl(MIN_REAL, MAX_REAL);
                    node = {std::to_string(val), nullptr, val, 0};
                }else if(rand_int <= 5){
                    node = {name_list[rand_int-2], operator_list[rand_int-2], 0.0, 2};
                }else if(rand_int <= NUM_OPERATORS + 1){
                    node = {name_list[rand_int-2], operator_list[rand_int-2], 0.0, 1};
                }

            }else{
                node = {"", nullptr, 0.0, -1};
            }

            nodes[i] = node;
            i++;

        }

    }

    return nodes;
}

vector<double> calcFunction(vector<Node> tree, vector<double>x_vec){

    int x_size = x_vec.size();
    int last_index = tree.size()-1;

    double a; double b; double y;
    vector<double> y_vec(x_size);

    for (int i=0; i<x_size; i++){
        if (tree[1].name=="x") {
                y = x_vec[i];
        }else{

            for (int j=last_index; j>0; j--){

                if(tree[j].type==1){

                    if(tree[2*j].name=="x"){
                        a = x_vec[i];
                    }else{
                        a = tree[2*j].value;
                    }

                    tree[j].value = tree[j].func(a, 0.0);

                }else if(tree[j].type==2){

                    if(tree[2*j].name=="x"){

                        a = x_vec[i];
                    }else{
                        a = tree[2*j].value;
                    }

                    if(tree[2*j+1].name=="x"){
                        b = x_vec[i];
                    }else{
                        b = tree[2*j+1].value;
                    }

                    tree[j].value =  tree[j].func(a, b);

                }
            }
            y = tree[1].value;
        }
        y_vec[i] = y;

    }

    return y_vec;
}

vector<double> calcFunction_old(vector<Node> tree_origin, vector<double>x_vec){

    int tree_size = tree_origin.size();
    static vector<Node> tree(tree_size);
    for (int i=0; i<tree_size; i++){
        Node node = tree_origin[i];
        tree[i] = node;
    }

    static int imax = tree_size-1;

    struct {
        double a; double b;

        double operator()(double x)
        {
            if (tree[1].name=="x") return x;

            for (int i=imax; i>0; i--){

                if(tree[i].type==1){

                    if(tree[2*i].name=="x"){
                        a = x;
                    }else{
                        a = tree[2*i].value;
                    }

                    tree[i].value = tree[i].func(a, 0.0);

                }else if(tree[i].type==2){

                    if(tree[2*i].name=="x"){
                        a = x;
                    }else{
                        a = tree[2*i].value;
                    }

                    if(tree[2*i+1].name=="x"){
                        b = x;
                    }else{
                        b = tree[2*i+1].value;
                    }

                    tree[i].value =  tree[i].func(a, b);

                }
            }
            return tree[1].value;
        }
    }func;

    int x_size = x_vec.size();
    vector<double> y_vec(x_size);

    for (int i=0; i<x_size; i++){
        y_vec[i] = func(x_vec[i]);
    }

    return y_vec;

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

        vector<double> y_vec_estimate = calcFunction(pop[i].nodes, x_vec);

        err = 0.0;
        for (int j=0; j<N_PTS; j++){
            err += (y_vec_estimate[j] - y_vec[j]) * (y_vec_estimate[j] - y_vec[j]);
        }
        err = sqrt(err);

        if (std::isnan(err)){
            pop[i].fitness = std::numeric_limits<double>::infinity();
        }else{
            pop[i].fitness = sqrt(err);
        }
        //showFunction(pop[i], true);
    }
}


void output(vector<Tree> pop, bool verbose){

    int imax = (int) (XMAX-XMIN)/DX;
    vector<double> x_vec(imax);
    vector<double> y_vec(imax);
    vector<double> x_vec_ref(N_PTS);
    vector<double> y_vec_ref(N_PTS);

    Tree best_ind = pop[0];
    double best_fitness = best_ind.fitness;
    double worst_fitness = pop[N_POP-1].fitness;
    double range = worst_fitness - best_fitness;
    int max_size=0;
    for (int i=0; i<N_POP; i++){
        int heap_size = pop[i].nodes.size();
        if (max_size < heap_size){max_size = heap_size;}
    }

    for (int i=0; i<=imax; i++){x_vec[i]=XMIN+i*DX;}
    y_vec = calcFunction(best_ind.nodes, x_vec);

    for (int i=0; i<N_PTS; i++){
        x_vec_ref[i] = point_table[i][0];
        y_vec_ref[i] = point_table[i][1];
    }

    ofs_history << g << ", " << best_fitness << endl;

    if ( (g + 1) % PATH_OUT_FREQ == 0){
        for (int i=0; i<imax; i++){
            ofs_bestpath << x_vec[i] << ", " << y_vec[i] << endl;
        }
    }

    if (g % DIST_OUT_FREQ == 0){
        for(int i=0; i<N_POP; i++){
            ofs_distribution << g << ", " << pop[i].fitness << endl;
        }
    }

    if (verbose){
        cout << "Generation: " << g
             << " / Fitness: " << best_fitness
             << " / Range: "   << range
             << " / max tree depth: " << log2(max_size)
             << endl;
    }


    if ( (g + 1) % PATH_OUT_FREQ == 0){
        string func_name = showFunction(best_ind, false);
        plt::clf();
        plt::xlim(0, 10);
        plt::ylim(-1, 1);
        plt::plot(x_vec_ref, y_vec_ref);
        plt::plot(x_vec, y_vec);
        plt::named_plot(func_name, x_vec, y_vec);
        plt::legend();
        plt::pause(0.1);
        //plt::show();
    }
}


int main()
{
    vector<Tree> pop;
    vector<Tree> elites(N_ELITES);
    vector<Tree> offspring;

    // Read data.
    std::string str;
    int i=0;
    while(getline(ifs, str)) {
		sscanf(str.data(), "%lf,%lf", &point_table[i][0], &point_table[i][1]);
        i++;
    }

    // Initial population
    pop = create_initial_pop(createTree);

    // Evaluation & sort.
    evaluate(pop);
    sort_pop(pop);

    // Output
    output(pop, true);

    for (g=1; g<NGEN; g++){

        // --- Selection
        elites = elitist_selection(pop);
        offspring = tournament_selection(pop);

        // --- Crossover
        crossover(oneptcx, offspring);

        // --- Mutation
        mutation(pm, offspring);

        // --- offspring + elites
        offspring.insert(offspring.end(), elites.begin(), elites.end());

         // --- Evaluation
        evaluate(offspring);
        sort_pop(offspring);

        // Update population
        pop = offspring;

        // --- Output
        output(pop, true);

    }

    return 0;

}

