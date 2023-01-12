// You may choose to use #include "CLI11.hpp" to pick up the repo version header file 
#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <string>

//
// Read in a copula file and re-format it with same convention as Formatter.
//
//                         Graham W. Wilson, 11-JAN-2023.
//


void Formatter(std::set<std::pair<std::pair<double,double>, int>> uvset,  std::set<std::pair<std::pair<double,double>, int>> vuset){

  // Process the sets and write out the empirical copula file 

  // The uvset and vuset have been sorted by either u or v given the ordering properties of std::set insertion.
  // Now make the empirical copula vector of pairs with the ranks where the events are already sorted by u using the 
  // event number as a key to get fast access to the corresponding v value.
  
    int N = uvset.size();
  
    std::vector<std::pair<int,int>> vfirst;          // For rank of u, event number
    vfirst.reserve(N);    
    std::map<int,int> mapSecond;                     // For event number (as key), rank of v

//    std::vector<std::pair<double,double>> vcopula;   // Empirical Copula with ranks re-scaled to (0,1)
//    vcopula.reserve(N);                              // May be useful to use vcopula directly, but prefer for now simply 
                                                       // writing this information to the copula file to decouple generation from fitting.   
      
    int rankuv=0;
    for (auto &uv : uvset) {
        rankuv++;
        vfirst.push_back(std::make_pair(rankuv, uv.second));
    }
    
    int rankvu=0;
    for (auto &vu : vuset) {
        rankvu++;
        mapSecond.insert(std::make_pair(vu.second, rankvu));        
    }

    std::ofstream fout;
    fout.open("CopulaGen.EDAT");
    fout.precision(10);   
         
    for (auto &u : vfirst){
        int ranku = u.first;
        int key = u.second;
        int rankv = mapSecond[key];
        // Calculate the empirical copula variables given by the rank in the range {1,2,....N}.
        // It is not 100% clear what is the correct procedure. 
        // Much of the literature uses rank/(N+1), but I agree with Joe in Dependence Modeling with Copulas (2014)  
        // that it is better to use (rank - 0.5)/N to populate more uniformly on the unit square. 
        // Both are motivated by mitigation of potential numerical issues near 0 or 1.
        std::pair<double, double> p = std::make_pair( (double(ranku)-0.5)/double(N), (double(rankv)-0.5)/double(N) );         
//        vcopula.push_back(p);
        fout << std::scientific << p.first << "  " << std::scientific << p.second << std::endl;        
    }
    
    fout.close();

}

void Reader(std::string copulafile="copula_4X5X-V3.EdatS"){

    std::set<std::pair<std::pair<double,double>, int>> uvset;
    std::set<std::pair<std::pair<double,double>, int>> vuset; 
    double u,v;
    
    std::cout << "Reading copula file named: " << copulafile << std::endl;

    std::ifstream copfile;
    copfile.open(copulafile);
    
    int i=0;
    while (copfile >> u >> v){
        i++;
        std::pair<double,double> puv= std::make_pair(u,v);
        std::pair<double,double> pvu= std::make_pair(v,u);
        uvset.insert(std::make_pair(puv,i));
        vuset.insert(std::make_pair(pvu,i));
    }
    std::cout << "Read " << i << " pairs" << std::endl;
    
    copfile.close();
    
    // Process the sets and write out the reformatted empirical copula file   
    Formatter(uvset, vuset);
    
}    

int main(int argc, char **argv) {

    CLI::App app{"Reformat emprical copula file"};
    
    int nevents = 9;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 13579L;
    app.add_option("-s,--seed", seed, "Seed");
           
    CLI11_PARSE(app, argc, argv);
    
    Reader();
       
    return 0;
    
}
