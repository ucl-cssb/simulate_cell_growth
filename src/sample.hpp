#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <map>

#include "cell.hpp"

using namespace std;


class Sample{
public:
    map<int, double> mut_vaf;
    double cellularity;
    double read_depth;
    double detect_limit;

    Sample();
    Sample(Clone clone, double cellularity, double read_depth, double detect_limit, mt19937 eng);
    ~Sample() = default;
    Sample(const Sample& other) = default;
    Sample(Sample&& other) = default;
    Sample& operator=(const Sample& other) = default;
    Sample& operator=(Sample&& other) = default;
    
    void print_vaf(string outfile);
};


#endif
