#include <vector>
#include <random>
#include <algorithm>
#include "utility.h"
using std::vector;

std::random_device seed_gen;
std::mt19937 mt_engine(seed_gen());
//std::mt19937 mt_engine(3);
int get_rand_range_int(int min_val, int max_val) {
    std::uniform_int_distribution<int> gen_rand_uni_int( min_val, max_val );
    return gen_rand_uni_int(mt_engine);
}

double get_rand_range_dbl(double min_val, double max_val) {
    std::uniform_real_distribution<double> gen_rand_uni_real( min_val, max_val );
    return gen_rand_uni_real(mt_engine);
}

int findIndex( vector<int> vec, int value ){
    vector<int>::iterator iter = std::find( vec.begin(), vec.end(), value);
    size_t index = std::distance( vec.begin(), iter );
    if(index == vec.size())
        {
            return -1;
        }
    return index;
}

