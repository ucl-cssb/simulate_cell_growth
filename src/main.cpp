#include <cstdlib>
#include <sstream>
#include <string>

#include "cell.hpp"
#include "sample.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
    Clone clone;

    int num_subclone = 0;
    int Nend = 1000;
    double birth_rate = log(2);
    double death_rate = 0;

    // Mutation rate per division per genome
    double mutation_rate = 0.01;

    string outdir = "./";
    double seed = 1;
    int verbose = 0;

    int num_clonal_mutation = 100;
    double min_clone_freq = 0.2;
    double max_clone_freq = 0.5;
    double tmin = 2;
    double tmax = 5;

    double read_depth =100;
    double cellularity = 0.6;
    double detect_limit = 0.05;


    if(argc > 1) num_subclone = atoi(argv[1]);
    if(argc > 2) Nend = atoi(argv[2]);
    if(argc > 3) birth_rate = atof(argv[3]);
    if(argc > 4) death_rate = atof(argv[4]);

    // Mutation rate per division per genome
    if(argc > 5) mutation_rate = atof(argv[5]);
    if(argc > 6) outdir = argv[6];
    if(argc > 7) seed = atof(argv[7]);
    if(argc > 8) verbose = atoi(argv[8]);

    if(argc > 9) num_clonal_mutation = atof(argv[9]);
    if(argc > 10) min_clone_freq = atof(argv[10]);
    if(argc > 11) max_clone_freq = atof(argv[11]);
    if(argc > 12) tmin = atof(argv[12]);
    if(argc > 13) tmax = atof(argv[13]);

    if(argc > 14) read_depth = atoi(argv[14]);
    if(argc > 15) cellularity = atoi(argv[15]);
    if(argc > 16) detect_limit = atoi(argv[16]);

    // net growth rate
    double lambda = birth_rate - death_rate;

    double tend = log(Nend)/(lambda); // The latest time that a subclone occurs

    vector<double> time_occur;
    vector<double> fitness;
    mt19937 eng;
    if(seed <= 0){
        eng = mt19937(time(0));
    }
    else{
        eng = mt19937(seed);
    }

    cout << "Net growth rate: " << lambda << endl;
    cout << "Simulation finish time (tumor doublings): " << tend << endl;

    if(tmax > tend){
        tmax = tend;
    }
    if(num_subclone > 0){
        // randomly select times that subclones occur
        // To simulate a subclone caused by WGD, use tend/2 since it often occurs early.
        uniform_real_distribution<double> runift(tmin, tmax);
        cout << "Subclone occurring time (tumor doublings): " << endl;
        for(int i=0; i<num_subclone; i++){
            time_occur.push_back(runift(eng));
            cout << "\t" << time_occur[i];
        }
        cout << endl;

        // randomly select fitness of subclones
        double s1, s2;
        if(min_clone_freq == 0){
            s1 = 0;
        }
        else{
            s1 = clone.get_subclone_fitness(lambda, min_clone_freq, tend, time_occur[0]);
        }
        if(max_clone_freq == 1){
            s2 = 1;
        }
        else{
            s2 = clone.get_subclone_fitness(lambda, max_clone_freq, tend, time_occur[0]);
        }
        cout << "Min and max subclone fitness: " << s1 << "\t" << s2 << endl;
        uniform_real_distribution<double> runiff(s1, s2);
        cout << "Subclone fitness and expected frequency: " << endl;
        double fexp = 0;
        for(int i=0; i<num_subclone; i++){
            fitness.push_back(runiff(eng));
            cout  << fitness[i] << ";";
            fexp = clone.get_subclone_freq_exp(lambda, fitness[i], tend, time_occur[0]);
            cout << fexp << "\t";
        }
        cout << endl;
    }


    cout << "Simulating tumor growth" << endl;
    cout << "Effective mutation rate (μ/β): " << mutation_rate/((birth_rate - death_rate)/birth_rate) << endl;

    clone.grow(num_subclone, num_clonal_mutation, fitness, time_occur, Nend, birth_rate, death_rate, mutation_rate, eng, verbose);

    string sep = "-";
    string ftype = ".txt";
    string midix = sep + to_string(num_subclone) + sep + to_string(Nend) + sep + to_string(int(mutation_rate)) + sep + to_string(int(seed));
    string suffix = midix + ftype;
    string outfile = "";

    if(verbose == 1){
        cout << "Printing out cells in a tumor clone" << endl;
        outfile = outdir + "all_cell_lineage"  + suffix;
        // cout << "file name" << outfile << endl;
        clone.print_lineage(clone.cells, outfile);

        outfile = outdir + "curr_vaf" + suffix;
        cout << "Computing VAF" << endl;
        map<int, double> mut_freq = clone.get_allele_freq();
        clone.print_map(mut_freq, outfile);
    }


    outfile = outdir + "subclone_freq" + suffix;
    cout << "Computing subclone frequency" << endl;
    map<int, double> subclone_freq = clone.get_subclone_freq();
    clone.print_map(subclone_freq, outfile);
    cout << clone.subclone_ID.size() << " subclone(s)" << endl;

    outfile = outdir + "subclone_lineage" + suffix;
    cout << "Printing subclone lineage" << endl;
    clone.print_clone_lineage(outfile);


    cout << "Sampling from the tumor" << endl;
    // double cellularity = 1;
    // double read_depth = 100;
    // // Assume detect_limit*read_depth reads are needed to call a variant
    // double detect_limit = 5 / read_depth;
    Sample sample(clone, cellularity, read_depth, detect_limit, eng);

    outfile = outdir + "sample_vaf" + suffix;
    sample.print_vaf(outfile);
    if(num_subclone>0){
        map<int, double> subclone_freq = clone.get_subclone_freq();
    }

    // Print details about the growing clone
    outfile = outdir + "summary" + suffix;
    cout << "Printing summary" << endl;
    clone.print_summary(outfile);

    return 0;
}
