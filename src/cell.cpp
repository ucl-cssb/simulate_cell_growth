#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
// #include <utility>
#include <random>
#include <cmath>
// #include <ctime>
#include <fstream>

#include "cell.hpp"

using namespace std;


Mutation::Mutation(){
        mut_ID = 0;
        time_occur = 0;
        vaf = 0;
        number = 1;
}


Mutation::Mutation(int mut_ID, double time_occur){
        this->mut_ID = mut_ID;
        this->time_occur = time_occur;
        this->vaf = 0;
        this->number = 1;
}


// Cell::Cell() {
//         cell_ID = 0;
//         parent_ID = 0;
//         clone_ID = 0;
//
//         birth_rate = log(2);
//         death_rate = 0;
//         fitness = 0;
//
//         mutation_rate = 0;
//
//         ploidy = 2;
//         num_division = 0;
//         time_occur = 0;
//         flag = 0;
// }


Cell::Cell(int cell_ID, int parent_ID, int clone_ID) {
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;
        this->clone_ID = clone_ID;

        this->birth_rate = log(2);
        this->death_rate = 0;
        this->fitness = 0;

        this->mutation_rate = 0;

        this->ploidy = 2;
        this->num_division = 0;
        this->time_occur = 0;
        this->flag = 0;
}


Cell::Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, double mutation_rate, double ploidy, double time_occur, double fitness = 0){
        this->cell_ID = cell_ID;
        this->parent_ID = parent_ID;

        this->birth_rate = birth_rate;
        this->death_rate = death_rate;
        this->fitness = fitness;

        this->mutation_rate = mutation_rate;

        this->ploidy = ploidy;
        this->time_occur = time_occur;
        this->flag = 0;
}


/*
mut_ID -- the ID of last mutation
*/
int Cell::generate_mutations(mt19937 eng, double mutation_rate, int& mut_ID, double time_occur){
        poisson_distribution<int> pois(mutation_rate);
        int nu = pois(eng);
        // cout << "Generating " << nu << " mutations" << endl;
        for (int j=0; j < nu; j++) {
                mut_ID += 1;
                // cout << mut_ID << "\t" << time_occur << endl;
                Mutation mut(mut_ID, time_occur);
                this->mutations.push_back(mut);
        }
        // cout << this->mutations.size() << endl;
        return nu;
}

/*
Generating specified number of mutations
*/
void Cell::generate_mutations_fixed(int& mut_ID, double time_occur, int num_mut){
        for (int j=0; j < num_mut; j++) {
                mut_ID += 1;
                // cout << mut_ID << "\t" << time_occur << endl;
                Mutation mut(mut_ID, time_occur);
                this->mutations.push_back(mut);
        }
        // cout << this->mutations.size() << endl;
}


// Some mutations may have multiple copies
int Cell::get_num_mut(){
    int sum = 0;
    for (auto mut : this->mutations){
        sum += mut.number;
    }
    return sum;
}


Cell* Cell::get_parent(vector<Cell>& cells){
    for(int i = 0; i < cells.size(); i++){
        Cell* cell = &cells[i];
        if(cell->cell_ID == parent_ID) return cell;
    }
    return NULL;
}


Clone::Clone(){
        clone_ID = 0;
        num_clonal_mutation = 0;
        time_occur = 0;
        frequency = 0;
        fitness = 0; // neutral evolution
}


/*
Initialize the start status of the clone
mut_ID: used to keep tracking the ID of new mutations, increased by 1 with a new mutation
nu: used to keep tracking the number of new mutations
*/
void Clone::initialize(double death_rate, double mutation_rate, int& mut_ID, int& nu, int num_clonal_mutation, int verbose = 0){
        this->clone_ID = 1;

        this->cells.clear();
        this->curr_cells.clear();

        this->subclone_ID.clear();
        this->subclone_parent.clear();
        this->subclone_time.clear();
        this->subclone_fitness.clear();
        this->num_clonal_mutation = num_clonal_mutation;

        // set parent as 0 to indicate MRCA
        Cell ncell(1, 0, 1);
        ncell.mutation_rate = mutation_rate;
        ncell.death_rate = death_rate;
        // Generate clonal mutations
        if(num_clonal_mutation > 0) {
            ncell.generate_mutations_fixed(mut_ID, 0, num_clonal_mutation);
            nu = num_clonal_mutation;
            // mut_ID += nu;
        }
        else if (mutation_rate > 0) {
            nu = ncell.generate_mutations(mt19937(mut_ID), mutation_rate, mut_ID, 0);
            // mut_ID += nu;
        }
        else{

        }

        if(verbose==1) {  // record all cells in the lineage
                this->cells.push_back(ncell);
        }
        this->curr_cells.push_back(ncell);
}


double Clone::get_rmax(){
    // Find  maximum birth rate (bmax) and maximum death rate (dmax) in the population
    double bmax = 0;
    double dmax = 0;
    for(unsigned int i=0; i<curr_cells.size(); i++) {
            Cell ci = curr_cells[i];
            double bi = ci.birth_rate;
            double di = ci.death_rate;
            if (bi > bmax) {
                    bmax = bi;
            }
            if(di > dmax) {
                    dmax = di;
            }
    }
    // cout << "Maximum birth rate: " << bmax << endl;
    // cout << "Maximum death rate: " << dmax << endl
    return bmax + dmax;
}

/*
   This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
   intput:
    Nend -- the number of cells in the final population;
    mutation_rate -- the Mutation rate per Cell division;
    time_occur -- defined in terms of population doublings
   output:
    a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
 */
void Clone::grow(int num_subclone, int num_clonal_mutation, vector<double> fitness, vector<double> time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, mt19937 eng, int verbose = 0){
    // time is defined in terms of population doublings
    double t = 0;
    int cell_count = 0;     // To count the total number of cells in histor
    // Initialize the simulation with one Cell
    int mut_ID = 0;
    int nu = 0; // The number of new mutations

    unsigned int fitmutant = 0;
    vector<double> timeN_occur;
    double lambda = birth_rate - death_rate;
    // convert time_occur from tumor doublings time to real time
    cout << "Subclone occurring time (tumor cell number):" << endl;
    for (unsigned int i = 0; i< time_occur.size(); i++){
        timeN_occur.push_back(ceil(exp(lambda * time_occur[i])));
        cout << "\t"  << timeN_occur[i];
    }
    cout << "\n";

    initialize(death_rate, mutation_rate, mut_ID, nu, num_clonal_mutation, verbose);
    cell_count += 1;
    while(this->curr_cells.size() < Nend) {
            if (this->curr_cells.size() == 0) {
                    t = 0;
                    mut_ID = 0;
                    nu = 0;
                    cell_count = 0;
                    initialize(death_rate, mutation_rate, mut_ID, nu, num_clonal_mutation, verbose);
                    cell_count += 1;
                    continue;
            }
            // Choose a random Cell from the current population
            uniform_int_distribution<> iunif(0, this->curr_cells.size() - 1);
            int rindex = iunif(eng);
            Cell rcell = this->curr_cells[rindex];
            // cout << "Selecting " << rindex+1 << "th cell" << endl;
            // cout << "Selecting cell " << rcell.cell_ID << endl;
            double rmax = get_rmax();

            // increase time
            uniform_real_distribution<double> runifu(0, 1);
            double tau = -log(runifu(eng)); // an exponentially distributed random variable
            double deltaT = tau/(rmax * this->curr_cells.size());
            t += deltaT;

            // draw a random number
            uniform_real_distribution<double> runif(0.0, rmax);
            double r = runif(eng);
            // cout << "random number " << r << endl;
            // birth event if r<birthrate
            if(r < rcell.birth_rate) {
                    // cout << "Number of generated cells " << this->cells.size() << endl;
                    // increase one Cell
                    // cout << "birth event at time " << t << endl;
                    int parent_ID = rcell.cell_ID;
                    // cout << "   parent " << parent_ID << endl;
                    Cell dcell1 = Cell(rcell);
                    dcell1.cell_ID = cell_count + 1;
                    dcell1.parent_ID = parent_ID;
                    dcell1.num_division = rcell.num_division + 1;

                    Cell dcell2 = Cell(rcell);
                    dcell2.cell_ID = cell_count + 2;
                    dcell2.parent_ID = parent_ID;
                    dcell2.num_division = rcell.num_division + 1;
                    cell_count += 2;
                    // cout << "   children " << dcell1.cell_ID  << "\t" << dcell2.cell_ID << endl;
                    // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                    if (mutation_rate>0) {
                        nu +=  dcell1.generate_mutations(mt19937(mut_ID), mutation_rate, mut_ID, t);
                        // mut_ID += nu;
                        nu +=  dcell2.generate_mutations(mt19937(mut_ID), mutation_rate, mut_ID, t);
                        // mut_ID += nu;
                    }
                    // introduce a fitter mutatant
                    if(fitmutant < num_subclone && this->curr_cells.size() >= timeN_occur[fitmutant]) {
                        double fitval = fitness[fitmutant];
                        cout << "Introducing a fitter mutatant with fitness " << fitval << endl;
                        dcell1.fitness = fitval;
                        dcell1.death_rate = runifu(eng) * rcell.death_rate;
                        dcell1.birth_rate = (1 + dcell1.fitness) * (rcell.birth_rate - rcell.death_rate) + dcell1.death_rate;
                        int clone_ID = this->clone_ID + fitmutant + 1;
                        cout << "   birth rate: " << dcell1.birth_rate << endl;
                        cout << "   death rate: " << dcell1.death_rate << endl;
                        dcell1.clone_ID = clone_ID;
                        this->subclone_ID.push_back(clone_ID);
                        this->subclone_parent[clone_ID] = rcell.clone_ID;
                        this->subclone_fitness[clone_ID] = fitval;
                        this->subclone_time[clone_ID] = t;
                        this->subclone_psize[clone_ID] = this->curr_cells.size();
                        fitmutant += 1;
                    }
                    if(verbose==1) {
                        for(unsigned int i = 0; i < this->cells.size(); i++) {
                                if (this->cells[i].cell_ID==rcell.cell_ID) {
                                        rcell.num_division += 1;
                                        rcell.flag = 1;
                                        this->cells[i] = rcell;
                                        break;
                                }
                        }
                    }
                    // Remove the parent cell from the list of current cells
                    // cout << "Removing cell " << this->curr_cells[rindex].cell_ID << endl;
                    this->curr_cells.erase(this->curr_cells.begin()+rindex);
                    dcell1.time_occur = t;
                    dcell2.time_occur = t;
                    this->curr_cells.push_back(dcell1);
                    this->curr_cells.push_back(dcell2);

                    if(verbose==1) {
                            this->cells.push_back(dcell1);
                            this->cells.push_back(dcell2);
                    }
            }
            // death event if b<r<b+d
            if(r >= rcell.birth_rate && r < rcell.birth_rate + rcell.death_rate) {
                    // cout << " death event" << endl;
                    if(verbose==1) {
                            for(unsigned int i = 0; i < this->cells.size(); i++) {
                                    if (this->cells[i].cell_ID==rcell.cell_ID) {
                                            rcell.flag = -1;
                                            this->cells[i] = rcell;
                                            break;
                                    }
                            }
                    }
                    this->curr_cells.erase(this->curr_cells.begin()+rindex);
            }
            // cout << "===========================" << endl;
    }
    // this->set_ploidy(this->curr_cells);
    // cout << "The average ploidy of this clone: " << this->ploidy << endl;
    cout << "Generated " << cell_count << " cells with " << nu << " mutations" << endl;
    cout << "End time: " << t << endl;
}


/*
   This method computes the average ploidy of a tumor population.
 */
double Clone::get_avg_ploidy(){
        double ploidy = 0;
        int num_cell = this->curr_cells.size();
        // Collect mutations
        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                ploidy += cell.ploidy;
        }
        ploidy = ploidy / num_cell;
        return ploidy;
}


/*
   This method computes the allele frequency of each mutation
 */
map<int, double> Clone::get_allele_freq(){
        map<int, double> mut_freq;

        int num_cell = this->curr_cells.size();
        // Collect mutations
        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                for(auto mut : cell.mutations) {
                        mut_freq[mut.mut_ID] += mut.number;
                }
        }
        // Compute frequency
        double ploidy = get_avg_ploidy();
        // cout << "Average ploidy of the population: " << ploidy << endl;
        for(auto it : mut_freq) {
                double vaf = it.second/num_cell;
                vaf = vaf / ploidy;
                mut_freq[it.first] = vaf;
        }
        return mut_freq;
}


/*
   This method computes the number of unique mutations in each subclone.
 */
map<int, int> Clone::get_subclone_nmut(){
        map<int, set<double>> subclone_muts;
        map<int, int> subclone_nmut;

        int num_cell = this->curr_cells.size();

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                for (auto mut : cell.mutations){
                    subclone_muts[cell.clone_ID].insert(mut.mut_ID);
                }
        }

        for(auto muts : subclone_muts) {
                subclone_nmut[muts.first] = muts.second.size();
        }

        return subclone_nmut;
}

/*
   This method computes the maximum number of cell division in each subclone.
 */
map<int, int> Clone::get_subclone_ndiv(){
        map<int, set<int>> subclone_divs;
        map<int, int> subclone_ndiv;

        int num_cell = this->curr_cells.size();

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                subclone_divs[cell.clone_ID].insert(cell.num_division);
        }

        // cout << "Number of unique divisions in subclones: " << endl;
        for(auto divs : subclone_divs) {
            // for(auto num : divs.second){
            //     cout << num << "\t";
            // }
            // cout << endl;
            // set<int>::iterator min = divs.second.begin();
            set<int>::reverse_iterator max = divs.second.rbegin();
            subclone_ndiv[divs.first] = *max;
            // cout << divs.first << "\t" << *min << "\t" << *max << "\n";
        }
        return subclone_ndiv;
}


/*
   This method computes the average number of cell division in each subclone.
 */
map<int, double> Clone::get_subclone_adiv(){
        map<int, vector<int>> subclone_divs;
        map<int, double> subclone_adiv;

        int num_cell = this->curr_cells.size();

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                // cout << cell.cell_ID << "\t" << cell.clone_ID << "\t" << cell.num_division << endl;
                subclone_divs[cell.clone_ID].push_back(cell.num_division);
        }

        // cout << "Number of divisions in subclones: " << endl;
        for(auto divs : subclone_divs) {
            int sum = 0;
            for (auto num : divs.second)
                sum += num;
            subclone_adiv[divs.first] = sum / (divs.second).size();
            // cout << "\t" << divs.first << "\t"<< sum << "\t" << (divs.second).size() << "\n";
        }

        return subclone_adiv;
}


/*
   This method computes the subclone frequencies of a tumor population.
 */
map<int, double> Clone::get_subclone_freq(){
        map<int, double> subclone_freq;
        int num_cell = this->curr_cells.size();

        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = this->curr_cells[i];
                subclone_freq[cell.clone_ID] += 1;
        }

        for(auto freq : subclone_freq) {
                subclone_freq[freq.first] = freq.second / num_cell;
        }

        return subclone_freq;
}


/*
   lambda: The net growth rate of the background host population
   freq: The frequency of a subclone
 */
double Clone::get_subclone_fitness(double lambda, double freq, double tend, double t1){
    double s = (lambda * t1 + log(freq / (1 - freq))) / (lambda * (tend - t1));
    return s;
}


/*
This function computes the theoretical subclone frequency.
   lambda: The net growth rate of the background host population
   freq: The frequency of a subclone
 */
double Clone::get_subclone_freq_exp(double lambda, double fitness, double tend, double t1){
    double numerator = exp(lambda * (1 + fitness) * (tend - t1));
    // double f = numerator / (numerator + exp(lambda * tend) + exp(lambda * (tend - t1)));
    double f = numerator / (numerator + exp(lambda * tend));
    return f;
}


void Clone::print_summary(string outfile) {
    ofstream out;
    out.setf(ios::fixed);
    out.setf(ios::showpoint);
    out.precision(9);
    out.open(outfile);

    map<int, int> subclone_nmut = get_subclone_nmut();
    map<int, double> subclone_freq = get_subclone_freq();
    map<int, int> subclone_ndiv = get_subclone_ndiv();
    map<int, double> subclone_adiv = get_subclone_adiv();

    double lambda = log(2);
    out << "Informaton for host population:" << endl;
    for(auto cell : curr_cells){
        if(cell.clone_ID == 0){
            lambda = cell.birth_rate - cell.death_rate;
            out << "\tMutation rate: " << cell.mutation_rate << endl;
            out << "\tBirth rate: " << cell.birth_rate << endl;
            out << "\tDeath rate: " << cell.death_rate << endl;
            out << "\tEffective mutation rate (μ/β): " << cell.mutation_rate / ((cell.birth_rate-cell.death_rate)/cell.birth_rate) << endl;
            out << "\tNumber of clonal mutation: "<< num_clonal_mutation << endl;
            break;
        }
    }

    int num_subclone = subclone_ID.size();
    out << "Number of subclones: " << num_subclone << endl;
    if (num_subclone > 0) {
            out << "Informaton for each subclone:" << endl;
            for(int i = 0; i < num_subclone; i++){
                out << "Subclone " << subclone_ID[i] << endl;
                for(auto cell : curr_cells){
                    if(cell.clone_ID == subclone_ID[i]){
                        out << "\tMutation rate: " << cell.mutation_rate << endl;
                        out << "\tBirth rate: " << cell.birth_rate << endl;
                        out << "\tDeath rate: " << cell.death_rate << endl;
                        out << "\tEffective mutation rate (μ/β): " << cell.mutation_rate / ((cell.birth_rate-cell.death_rate)/cell.birth_rate) << endl;
                        break;
                    }
                }
                out << "\tFrequency: " << subclone_freq[subclone_ID[i]] << endl;
                out << "\tNumber of mutations in subclone: " << subclone_nmut[subclone_ID[i]] << endl;
                out << "\tFitness advantage: " << subclone_fitness[subclone_ID[i]] << endl;
                out << "\tTime subclone emerges (simulation time): " << subclone_time[subclone_ID[i]] << endl;
                out << "\tNumber of divisions: " << subclone_ndiv[subclone_ID[i]] << endl;
                out << "\tAverage number of divisions per cell: " << subclone_adiv[subclone_ID[i]] << endl;
                out << "\tPopulation size when subclone emerges: " << subclone_psize[subclone_ID[i]] << endl;
                out << "\tTime subclone emerges (tumor doublings): " << log(subclone_psize[subclone_ID[i]])/(lambda) << endl;
                out << "\tParent of subclone (0 is host): " << subclone_parent[subclone_ID[i]];
                out << endl;
            }
    }
    else{
            out << "No clones, tumour growth was neutral\n\n";
    }
    out.close();
}


/*
   This method prints out the details of cells in a clone
 */
void Clone::print_lineage(vector<Cell> cells, string outfile){
        ofstream out;
        out.open(outfile);
        int num_cell = cells.size();
        cout << "There are " << num_cell << " cells in the history of tumor growth" << endl;
        // string header = "id\tparent_ID\tflag\tbirth_rate\tdeath_rate\tmutation_rate\tploidy\tnum_division\ttime_occur\n";
        string header = "id\tparent_ID\tflag\tnum_mut\tclone_ID\ttime_occur\tbranch_len\tfitness\n";
        out << header;
        for(unsigned int i = 0; i < num_cell; i++) {
                Cell cell = cells[i];
                // cout << cell.cell_ID  << "\t" << cell.parent_ID << "\t" << cell.flag << "\t" << cell.birth_rate << "\t" << cell.death_rate << "\t" << cell.mutation_rate  << "\t" << cell.ploidy  << "\t" << cell.num_division << "\t" << cell.time_occur << endl;
                int num_mut = cell.get_num_mut();
                // int num_mut = cell.mutations.size();
                if(cell.parent_ID == 0) continue;
                Cell* pcell = cell.get_parent(cells);
                double branch_len = cell.time_occur - pcell->time_occur;
                double fitness = 1;
                if(cell.clone_ID > 1){
                  fitness = subclone_fitness[cell.clone_ID];
                }
                // write the cell lineages in a file for visualization
                out << cell.cell_ID << "\t" << cell.parent_ID << "\t" << cell.flag << "\t" << num_mut << "\t" << cell.clone_ID << "\t" << cell.time_occur << "\t" << branch_len << "\t" << fitness << endl;
                // cout << "\t" << num_mut << " mutations in the cell: " << endl;
                // for(unsigned int j = 0; j < num_mut; j++){
                //     Mutation mut = cell.mutations[j];
                //     cout << "\t" << mut.mut_ID << "\t" << mut.time_occur << endl;
                // }
        }
        out.close();
}


/*
   This method prints out the relationship among clones
 */
void Clone::print_clone_lineage(string outfile){
        cout << "There are " << this->subclone_ID.size() << " clones in the history of tumor growth" << endl;
        ofstream out;
        out.open(outfile);
        string header = "ID\tparent_ID\n";
        out << header;
        for (auto id : this->subclone_ID) {
          out << id << "\t" << this->subclone_parent[id] << endl;
        }
        out.close();
}
