#!/bin/sh

WEIGHT=$1

./BetaModeler -n 1000000 -w ${WEIGHT} --par1 12.61 --par2 -0.64205 --par3 38.9 --par4 -0.3891 -f 0.939

exit
