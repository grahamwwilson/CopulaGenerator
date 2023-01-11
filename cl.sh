#!/bin/sh
#
# You likely will need to adapt the CLI11.hpp include (it is provided) 
# or install the CLI11 header file on your system. 
# CLI11 provides a light-weight C++ method for command line argument parsing.
#
# Once you have that working, simply 
# > ./cl.sh 
# should be sufficient to make the Generate executable as this is the default.
#

fn=${1:-Generate}

# Compile with high optimization. No need for any libraries with this for now.
g++ -O3 ${fn}.cpp -o ${fn}

exit
