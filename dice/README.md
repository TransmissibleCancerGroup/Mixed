# Dice - repeat mutation simulator

Usage: dice --threads 2 --reps 1000 --mutations 15000 spectrum.txt opps.txt

The above command runs 1000 simulation repetitions.
Each repetition generates 15000 mutations.
Each mutation is:
  - Randomly given a type according to probabilities in `spectrum.txt`
  - Placed uniformly at random in a position according to mutational
    opportunities in `opps.txt`

Output: the proportion of repetitions that generated a repeated mutation
(same mutation, same position)

The output is an estimate of the probability that there will be at least
one recurring mutation among a set of 15000.


# Compiling the program

dice is written in D. It can be compiled using dub:

    dub build

A fast version can be built with

    dub build --compiler=ldc2 --build=release

dice_cpp is written in C++. It can be build using cmake:

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    cp dice ../dice_cpp

