#!/bin/sh

N=$1

# Generate independent distributions (faster than option -c 2) and does 
# not depend on numerical inversion algorithm
time ./FullMonty -n ${N} -s 10011 -i >FM0.out
mv FullMonty.root FullMonty-${N}-C0.root

# Standard case with fitted copula model
time ./FullMonty -n ${N} -s 10001 >FM1.out
mv FullMonty.root FullMonty-${N}-C1.root

# Use independence copula
time ./FullMonty -n ${N} -s 10021 -c 2 >FM2.out
mv FullMonty.root FullMonty-${N}-C2.root

exit
