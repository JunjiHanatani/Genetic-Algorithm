#include "Vector3D.h"
#include <vector>
#include <algorithm>
//----------------------------------------------------
// Vector operation
//----------------------------------------------------

double calcDistance(vector<double>&p0, vector<double>&p1){
    vector<double> vec = sub(p0, p1);
    return calcNorm(vec);
}

double calcNorm(vector<double>&vec){
    double sum_of_squares = std::accumulate(vec.begin(), vec.end(), 0.0,
                   [](double sum, double x){return sum + x*x;});
    return sqrt(sum_of_squares);
}

vector<double> calcUnitVector(vector<double>&vec){
    double abs = 1.0/calcNorm(vec);
    return scaling(vec, abs);
}

vector<double> calcUnitVector2(vector<double>&p1, vector<double>&p2){
    vector<double> vec = sub(p2, p1);
    double abs = 1.0/calcNorm(vec);
    return scaling(vec, abs);
}

vector<double> calcVector(vector<double>&vec, double value){
    vector<double> unit_vec = calcUnitVector(vec);
    return scaling(unit_vec, value);
}

double dotProduct(vector<double>&vec1, vector<double>&vec2){
    vector<double> vec = mul(vec1, vec2);
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return sum;
}

vector<double> add(vector<double>& v1, vector<double>& v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x+y;});
    return v3;
};

vector<double> sub(vector<double>& v1, vector<double>& v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x-y;});
    return v3;
};

vector<double> mul(vector<double>& v1, vector<double>& v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x*y;});
    return v3;
};

vector<double> div(vector<double>& v1, vector<double>& v2){
    vector<double> v3(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                   [](double x, double y){return x/y;});
    return v3;
};

vector<double> scaling(vector<double>& vec, double value){
    vector<double> vec_res(vec.size());
    std::transform(vec.begin(), vec.end(), vec_res.begin(),
                   [&value](double x){return x*value;});
    return vec_res;
};
