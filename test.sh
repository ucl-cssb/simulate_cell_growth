#!/bin/bash

outdir=./test/
if [[ ! -d $outdir ]]; then
    mkdir -p $outdir
fi

num_subclone=0
num_cell=1000

birth_rate=0.6931472
death_rate=0
mutation_rate=5

seed=1
verbose=0

# Parameters related to subclone growth
# num_clonal_mutation=200
# min_clone_freq=0.5    # The minimal frequency of a subclone
# max_clone_freq=0.6
# tmin=2        # The earliest time that a subclone occurs
# tmax=3

# Parameters related to SNV sampling
# read_depth=200
# cellularity=1
# detect_limit=0.05


/usr/bin/time ./bin/sim $num_subclone $num_cell $birth_rate $death_rate $mutation_rate $outdir $seed $verbose
