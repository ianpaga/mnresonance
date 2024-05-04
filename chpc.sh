#!/bin/bash

# compile mu100-400by400by100_E***.cpp in the chosen start_$2 directory

cd start_$2/
c++ -O4 -I/groups/astro/ianpaga/boost_1_71_0 -fopenmp $1.cpp -o $1.out
##c++ -O4 -Wall -I/groups/astro/ianpaga/boost_1_71_0 -fopenmp $1.cpp -o $1.out
