#!/bin/sh
#
# This needs the CLI11.hpp include (it is provided) or you can install the CLI11 
# header file on your system. 
# CLI11 provides a light-weight C++ method for command line argument parsing.
#
# Simply 
# ./clroot.sh is sufficient to make the BivariateModeler executable as this is 
# the default.
#

fn=${1:-BivariateModeler}

# Compile with high optimization and link against ROOT libraries
g++ -O3 ${fn}.cpp -o ${fn} `root-config --cflags --glibs`

exit
