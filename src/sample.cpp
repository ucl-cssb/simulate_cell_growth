#include <iostream>
#include <fstream>
#include <map>
#include <random>

#include "sample.hpp"

using namespace std;

Sample::Sample(){
    cellularity = 1;
    read_depth = 100;
    detect_limit = 0.1;
}


Sample::Sample(Clone clone, double cellularity, double read_depth, double detect_limit, mt19937 eng){
    binomial_distribution<int> binom;
    poisson_distribution<int> pois(read_depth);

    int num_cell = clone.curr_cells.size();
    double seq_prob = read_depth/num_cell;
    map<int, double> mut_freq = clone.get_allele_freq();

    for(auto it : mut_freq){
        // ofile << it->first  << "\t" << it->second << endl;
        // binom = binomial_distribution<int>(num_cell, seq_prob);
        // double seq_depth = binom(eng);
        double seq_depth = pois(eng);
        // cout << "depth for mut " << it.first << " is " << seq_depth << endl;

        double vaf_true = it.second;
        double read_prob = vaf_true * cellularity;
        if(read_prob<detect_limit) continue;
        // cout << "    prob is " << read_prob << endl;
        binom = binomial_distribution<int>(seq_depth, read_prob);
        int read_count = binom(eng);

        double vaf = read_count/seq_depth;
        this->mut_vaf[it.first] = vaf;
    }
}


void Sample::print_vaf(string outfile){
    ofstream ofile;
    ofile.setf(ios::fixed);
    ofile.setf(ios::showpoint);
    ofile.precision(9);
    ofile.open(outfile);
    for(auto it=this->mut_vaf.begin(); it!=this->mut_vaf.end(); it++){
        ofile << it->first << "\t" << it->second << endl;
    }
    ofile.close();
}
