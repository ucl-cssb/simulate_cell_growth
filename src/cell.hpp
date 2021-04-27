#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>


using namespace std;

class Mutation
{
public:
    int mut_ID;
    int size;   // 1 for SNV
    int type;   // mutation type. 1 -- SNV, 2 -- CNV
    int start;
    int end;
    double time_occur;
    double vaf;
    int number;

    Mutation();
    Mutation(int mut_ID, double time_occur);
    ~Mutation() = default;
    Mutation(const Mutation& other) = default;
    Mutation(Mutation&& other) = default;
    Mutation& operator=(const Mutation& other) = default;
    Mutation& operator=(Mutation&& other) = default;
};


class Cell
{
public:
    int cell_ID;
    int parent_ID;
    int clone_ID;

    double birth_rate;
    double death_rate;
    double fitness;

    double mutation_rate;

    double ploidy;
    int num_division;
    double time_occur;
    int flag;   // whether the cell is alive or not. 0:new, 1:divided, -1: death

    vector<Mutation> mutations;

    // Cell();
    Cell(int cell_ID, int parent_ID, int clone_ID);
    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, double mutation_rate, double ploidy, double time_occur, double fitness);

    ~Cell() = default;
    Cell(const Cell& other) = default;
    Cell(Cell&& other) = default;
    Cell& operator=(const Cell& other) = default;
    Cell& operator=(Cell&& other) = default;

    int generate_mutations(mt19937 eng, double mu, int& mut_ID, double time_occur);
    void generate_mutations_fixed(int& mut_ID, double time_occur, int num_mut);

    void update_mut_count(int multiple);    // Increase the number of mutations due to WGD
    int get_num_mut();

    Cell* get_parent(vector<Cell>& cells);
};



/*
to represent a population of cells
*/
class Clone
{
public:
    int clone_ID;
    int num_clonal_mutation;
    vector<Cell>  cells;    // all the cells in the history, used for checking lineage history
    vector<Cell>  curr_cells;   // only available cells at present

    // variables to store subclone informaton
    vector<int> subclone_ID;  // ID of subclones
    map<int, double> subclone_time;   // Time subclone emerges
    map<int, double> subclone_fitness;
    map<int, int> subclone_parent;  // parents of subclones to track parent-child relationship of subclones
    map<int, int> subclone_psize;   // size of subclones

    double time_occur;
    double frequency;
    double fitness;     // fitness of entire clone population

    Clone();
    Clone(int clone_ID, double time_occur, double fitness);

    ~Clone() = default;
    Clone(const Clone& other) = default;
    Clone& operator=(const Clone& other) = default;
    Clone& operator=(Clone&& other) = default;

    // Methods related to cell growth in a clone
    void initialize(double death_rate, double mutation_rate, int& mut_ID, int& nu, int num_clonal_mutation, int verbose);
    double get_rmax();
    void grow(int num_subclone, int num_clonal_mutation, vector<double> fitness, vector<double> time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, mt19937 eng, int verbose);

    double get_avg_ploidy();
    map<int, double> get_allele_freq();

    // methods related to subclone growth
    double get_subclone_fitness(double lambda, double freq, double tend, double t1);
    double get_subclone_freq_exp(double lambda, double fitness, double tend, double t1);
    map<int, double> get_subclone_freq();
    map<int, int> get_subclone_nmut();
    map<int, int> get_subclone_ndiv();
    map<int, double> get_subclone_adiv();

    void print_lineage(vector<Cell> cells, string outfile);
    void print_clone_lineage(string outfile);
    void print_summary(string outfile);

    template <typename T>
    void print_map(map<int, T> m, string outfile){
        ofstream out;
        out.setf(ios::fixed);
        out.setf(ios::showpoint);
        out.precision(9);
        out.open(outfile);

        for(auto it : m) {
                out << it.first << "\t" << it.second << endl;
        }

        out.close();
    }
};

#endif
