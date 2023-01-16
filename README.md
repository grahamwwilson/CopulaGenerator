# CopulaGenerator
This packaged started life for the singular purpose of generating empirical 
copulas from a specified distribution. It now has some other components related 
to modeling 3-body decays.

## Executables
### Generate           
Empirical copula distribution generator (for specific model)
### Reformat           
Take (x,y) dataset and reformat as ordered empirical distribution based on ranks.
### BivariateModeler   
Sample (x,y) and related values from Nojiri, Yamada expressions for neutralino 3-body decay, A -> B f fbar 
### UnivariateModeler  
Sample neutralino energy for A -> B f fbar using Baer, Tata expression.

## Compilation
./cl.sh for compilation of Generate

./cl.sh Reformat for compilation of Reformat

./clroot.sh for compilation of BivariateModeler

./clroot.sh UnivariateModeler for compilation of UnivariateModeler


## Help
For help on run time parameters:

./Generate --help

./Reformat --help

./BivariateModeler --help

./UnivariateModeler --help
