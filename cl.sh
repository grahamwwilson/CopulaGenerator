#!/bin/sh
#
# You likely will need to change the CLI11.hpp include.
#

fn=${1:-Generate}

# Compile with high optimization. No need for any libraries with this for now.
g++ -O3 ${fn}.cpp -o ${fn}

exit
