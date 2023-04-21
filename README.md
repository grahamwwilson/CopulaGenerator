# CopulaGenerator
This package started life with the singular purpose of generating empirical 
copulas from a specified distribution. It now includes a related reformatting 
utility, and some other components related to sampling from 3-body decays.

## Executables
### Generate           
Empirical copula distribution generator (for specific model)
### Reformat           
Take (x,y) dataset and reformat as ordered empirical distribution based on ranks.
### BivariateModeler   
Sample (x,y) and related values from Nojiri, Yamada expressions for neutralino 3-body decay, A -> B f fbar 
### UnivariateModeler  
Sample neutralino energy for A -> B f fbar using Baer, Tata expression.
### FullMonty
Generate (x1, x2) from luminosity spectrum using 4-region approach 
(peak, arm1, arm2, body) with separate double beta distributions for 
the body and arms. Implementation includes Gaussian beam energy spread 
and a parametrization of the (x1,x2) dependence copula distribution for 
body events. In total, there are 3 region parameters, 2 resolution parameters, 
3 dependence parameters, and 5 shape parameters for the marginal double-beta 
distributions in the body and arms prior to beam energy spread convolution, 
for a total of 18 parameters.

## Compilation
./cl.sh for compilation of Generate

./cl.sh Reformat for compilation of Reformat

./clroot.sh for compilation of BivariateModeler

./clroot.sh UnivariateModeler for compilation of UnivariateModeler

./clroot.sh FullMonty for compilation of FullMonty

## Help
For help on run time parameters:

./Generate --help

./Reformat --help

./BivariateModeler --help

./UnivariateModeler --help

./FullMonty --help
