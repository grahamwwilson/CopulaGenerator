// You may choose to use #include "CLI11.hpp" to pick up the repo version header file 
#include "CLI11.hpp"
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <map>
#include <string>

typedef std::mt19937 RandomNumberGenerator;

std::string copulafile;

//
// Generate empirical copula using order statistics for Goodness-Of-Fit tests.
// The null hypothesis we wish to test is that the (E1,E2) dependence structure of the 
// GP4X5X data with deltax_min=2.0e-6 belongs to the Clayton + Ali-Mikhail-Haq mixture 
// copula distribution.
//
// C(u,v) = w C1(u,v) + (1-w) C2(u,v)
//
// where C1(u,v) = Clayton(thetaC) = pow ( pow(u,-thetaC) + pow(v,-thetaC) - 1, -1/thetaC ),
//       C2(u,v) = AMH(thA) = u*v/(1 - thA*(1-u)*(1-v) ), and 
//       w is the weight of the Clayton component.
//
// We will use the so-called parametric bootstrap procedure following 
// Appendix A of Genest, Remillard, Beaudoin, 
// Insurance Mathematics and Economics 44 (2009) 199-213. 
//
// We take the 3 parameters of this copula model using the estimated parameters 
// from the original pseudo-likelihood fit to data based on ranked observations to the 583,584 events, 
// and will generate multiple random samples according to the estimated copula parameters.
// Each random sample will then be fitted with the same model (with the 3 parameters floating) 
// in the same manner, and used to determine empirically the p-value based on 
// various statistics including Cramer-von-Mises, Kolomogorov, and Chi-Square whether 
// the fit model is adequate. 
// Note that given that we are fitting order-statistics based copulas where the marginals are 
// exactly uniform distributions, the reduction in degrees-of-freedom for for example 
// a chi-squared based GoF estimate is more than just the number of free parameters.
//
// Purpose of this code. 
// Generate a random sample file sorted as an empirical copula.
//
// Method.
// Fortunately both the Clayton and AMH copulas are relatively 
// straightforward to simulate using the distribution inverse method and the mixture can 
// be assured using Bernouilli trials.
// Note that as both are Archimedean copulas it is also possible to  
// use the generator procedure of Genest and MacKay (1986). I considered this  
// but concluded that this was not necessary given that the distribution inverse 
// method was tractable enough. It could be faster, but I doubt this is a bottle-neck.
// Additionally this code does the sorting and preparation of the empirical copula 
// using STL vectors, sets, and maps.
//
//                         Graham W. Wilson, 10-JAN-2023.
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
    fout.open(copulafile);
    std::cout << "Writing copula file named " << copulafile << std::endl;
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

void Generator(int nevents, unsigned long int seed){

// Hard code the best fit parameters based on fit to GP4X5X data (583,584 events).
// Fit parameters 
// Clayton dependence parameter 0.04969906094
// AMH dependence parameter     0.3598358846
// Clayton phiwt parameter      0.4436183234 wt = 0.8157784722

    double thetaC  = 0.04969906094;
    double thA     = 0.3598358846;
    double weightC = 0.8157784722; 

    RandomNumberGenerator g(seed);
    std::uniform_real_distribution<double> uniform; 
    
    std::set<std::pair<std::pair<double,double>, int>> uvset;
    std::set<std::pair<std::pair<double,double>, int>> vuset; 
    
    double u1,u2;
    
    for (int i=0; i<nevents; i++){
        double which = uniform(g); // Use Bernouilli trials with p=weightC to decide whether Clayton or AMH
        u1 = uniform(g);           // We always use this one as the first copula random variate.
        double v2 = uniform(g);    // v2 is used as an auxiliary random variable to solve for u2 using v2 = dC(u,v)/du = f(u=u1, v=u2)
        int idistbn;
        if ( which <= weightC ){
  // Clayton (also known by other names. Joe (2014) calls it the bivariate Mardia-Takahasi-Clayton-Cook-Johnson copula
            idistbn = 0;
            u2 = pow(1.0 - pow(u1, -thetaC) + pow( 1.0/(v2*pow(u1, 1.0+thetaC)), thetaC/(1.0+thetaC) ) , -1.0/thetaC);      
        }
        else{
  // AMH
            idistbn = 1;
            double A = v2*pow(thA*(1.0-u1),2.0) - thA;
            double B = 2.0*v2*thA*(1.0-u1)*(1.0 - thA*(1.0-u1) ) + thA - 1.0;
            double C = v2*(1.0 - 2.0*thA*(1.0-u1) + pow(thA*(1.0-u1),2.0)); 
            double D = B*B - 4.0*A*C;
            u2 = v2;
            if (D>=0.0){
                double xplus  = (-B + sqrt(D))/(2.0*A);
                double xminus = (-B - sqrt(D))/(2.0*A);
                u2 = xminus;
                if(u2<=0.0 || u2>=1.0){
                   std::cout << "Scream 1: A, B, C, D " << A << " " << B << " " << C << " " << D << " " << xplus << " " << xminus << std::endl;                   
                }
            }
            else{
               std::cout << "Scream 2: A, B, C, D " << A << " " << B << " " << C << " " << D << std::endl;
            } 
        }
//        fout << std::setw(10) << i+1 << "  " << std::scientific << u1 << "  " << std::scientific << u2 << "  " << idistbn << std::endl;

 // At this stage (u1,u2) are random variables on (0,1) from the copula distribution

        std::pair<double,double> puv=std::make_pair(u1,u2);
        uvset.insert(std::make_pair(puv,i+1));              // Will be sorted by increasing values of u1
        std::pair<double,double> pvu=std::make_pair(u2,u1);
        vuset.insert(std::make_pair(pvu,i+1));              // Will be sorted by increasing values of u2
        
    }
    
// Process the sets and write out the empirical copula file   
    Formatter(uvset, vuset);
    
}

int main(int argc, char **argv) {

    CLI::App app{"Empirical copula generator"};
    
    int nevents = 9;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 13579L;
    app.add_option("-s,--seed", seed, "Seed");

    std::string filename = "CopulaGen.EDAT";
    app.add_option("-o,--outputfile", filename, "Output copula file");

    CLI11_PARSE(app, argc, argv);

    std::cout << "nevents  " << nevents << std::endl;
    std::cout << "seed     " << seed << std::endl;
    std::cout << "filename " << filename << std::endl;

    copulafile = filename;             // Write filename to globally declared copulafile string
    
    Generator(nevents, seed);
       
    return 0;
    
}
