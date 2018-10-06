#ifndef _UTILITY_H_
#define _UTILITY_H_

// Random Seed
extern std::mt19937 mt_engine;

int get_rand_range_int(int, int);
double get_rand_range_dbl(double, double);

/*
class CSample
{
public:
    void set(int num);
    int get();
private:
    int m_num;
};
*/

#endif // _UTILITY_H_
