#include "CLI11.hpp"     // Here use #include "CLI11.hpp" to pick up the local header file (distributed with this repository) 
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11. This works if installed on system.
#include <cmath>
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <map>
#include <string>
typedef std::mt19937 RandomNumberGenerator;

//
// Generate 2-mixture beta distribution using hit and miss MC
//
//                         Graham W. Wilson, 16-JAN-2023.
//

// For inverting the cdf I think I'll need to do some numerical integration modeled on that of EnergyFits/RooMyConvolvedBetaPdf.cxx


double functionSecondDerivative(double t, double normFactor, double eta, double a2, double a3, double a2p, double a3p, double frac){
// 
// Evaluate the second derivative of the function at t for potential use in Halley's method.
//
// FIXME Still need to write this. See p21 of new logbook for derivation 17-JAN-2023.
    double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
    double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
    double tintegrand = frac*beta1 + (1.0-frac)*beta2;
    double pdf = tintegrand/normFactor;          

    return pdf;
    
}


double functionGradient(double t, double normFactor, double eta, double a2, double a3, double a2p, double a3p, double frac){
// 
// Evaluate the function gradient at t
//

    double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
    double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
    double tintegrand = frac*beta1 + (1.0-frac)*beta2;
    double pdf = tintegrand/normFactor;          

    return pdf;
    
}

std::pair<double, double> Evaluate(double tmin, double tmax, int Nbase, double eta, double a2, double a3, double a2p, double a3p, double frac){
// 
// Evaluate the un-normalized cumulative distribution function using Boole's rule.
// This is a 5-point closed Newton-Cotes formula. 
// (the first call will be to calculate the integral over the full range)
//
    double t;
    double integral = 0.0;
    double integrals = 0.0;  // Keep track of Simpson's rule estimate too for error estimation
    int N = 4*Nbase;         // Needs to be a multiple of 4
    double dt = (tmax-tmin)/double(N);
    
    for (int i = 0; i <=N; i++){
        t = tmin + double(i)*dt;
        double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
        double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
        double tintegrand = frac*beta1 + (1.0-frac)*beta2;          
        if ( i==0 || i==N){
            integral += 7.0*tintegrand;
            integrals += tintegrand;
        }
        else if (i%2==1){            
            integral += 32.0*tintegrand;
            integrals += 4.0*tintegrand;
        }    
        else if (i%4==2){
            integral += 12.0*tintegrand;
            integrals += 2.0*tintegrand;
        }
        else{
            integral += 14.0*tintegrand; // interior every other even abscissa (counting from 0: So i=4,8,12,..))     
            integrals += 2.0*tintegrand;
        }
    }
    integrals *= dt/3.0;
    integral *= 2.0*dt/45.0;
    std::pair<double, double> p = std::make_pair(integral, integrals);
    return p;
    
}

double rootFinderAlt( double tlo, double thi, double flo, double fhi, double u, double normFactor, double tmin, double tmax, 
                      double eta, double a2, double a3, double a2p, double a3p, double frac){
    //
    // Linear interpolation.  Two points (tlo, flo), (thi, fhi).
    // Line through the two points is f = gradient*t + intercept
    // Whose zero is at t = -gradient/intercept
    //
    // This estimate is achieved without any need for further function evaluation, 
    // but it does not carry with it error estimation.
    // Although maybe it is reasonable to assign an error related to the difference between the mean 
    // and the linear interpolation estimate.
    //
    double gradient = (fhi-flo)/(thi-tlo);
    double intercept = flo - gradient*tlo;
    double troot = -intercept/gradient;

    double tmid = (tlo+thi)/2.0;

    // Now try Newton-Raphson
    int niter;
    double tol = 1.0e-12;
    int Nsteps = 40;
    double error = 1.0;
    
    niter = -1;
    
    std::pair<double, double> dintegral;
    double cdfestimate;   
    
    while (niter < 5 && error > tol){
    
  // Only calculate the corrections to the points in the look up table (so can afford many fewer steps).
        if(troot > tmid){
            dintegral = Evaluate( troot, thi, Nsteps, eta, a2, a3, a2p, a3p, frac );
            cdfestimate = (fhi + u) - (dintegral.first/normFactor);
        }
        else{
            dintegral = Evaluate( tlo, troot, Nsteps, eta, a2, a3, a2p, a3p, frac );
            cdfestimate = (flo + u) + (dintegral.first/normFactor);            
        }    
            
        double fn = cdfestimate - u;
        double fngrad = functionGradient( troot, normFactor, eta, a2, a3, a2p, a3p, frac);
    // Update the root anyway at this point
        double dtroot = -(fn/fngrad);
        troot+= dtroot;
        error = std::abs(fn);
        niter++;
        std::cout << "NR iteration " << niter << " dtroot " << dtroot << " error " << error << std::endl;
 
    }    

    std::cout << "Returning after " << niter << " iterations " << std::endl;   
    
    return troot;
    
}

double rootFinder( double tlo, double thi, double flo, double fhi, double u, double normFactor, double tmin, double tmax, 
                   double eta, double a2, double a3, double a2p, double a3p, double frac){
    //
    // Linear interpolation.  Two points (tlo, flo), (thi, fhi).
    // Line through the two points is f = gradient*t + intercept
    // Whose zero is at t = -gradient/intercept
    //
    // This estimate is achieved without any need for further function evaluation, 
    // but it does not carry with it error estimation.
    // Although maybe it is reasonable to assign an error related to the difference between the mean 
    // and the linear interpolation estimate.
    //
    double gradient = (fhi-flo)/(thi-tlo);
    double intercept = flo - gradient*tlo;
    double troot;
    
    troot = -intercept/gradient;

    // Now try Newton-Raphson
    int niter;
    double tol = 1.0e-12;
    int Nsteps = 200;
    double error = 1.0;
    
    niter = -1;
    
    std::pair<double, double> integral;
    double cdfestimate;   
    
    while (niter < 5 && error > tol){
    
        if(u < 0.5){
            integral = Evaluate( tmin, troot, Nsteps, eta, a2, a3, a2p, a3p, frac );
            cdfestimate = integral.first/normFactor;
        }
        else{
            integral = Evaluate( troot, tmax, Nsteps, eta, a2, a3, a2p, a3p, frac );
            cdfestimate = 1.0 - (integral.first/normFactor);            
        }    
        double fn = cdfestimate - u;
        double fngrad = functionGradient( troot, normFactor, eta, a2, a3, a2p, a3p, frac);
    // Update the root anyway at this point
        double dtroot = -(fn/fngrad);
        troot+= dtroot;
        error = std::abs(fn);
        niter++;
        std::cout << "NR iteration " << niter << " dtroot " << dtroot << " error " << error << std::endl;
 
    }    

    std::cout << "Returning after " << niter << " iterations " << std::endl;   
    
    return troot;
    
}

void UniGenerator(int nevents, unsigned long int seed, double dxmin, double eta, double a2, double a3, double a2p, double a3p, double frac, double WTMAX){

    double maxweight;
    
    double tmin = std::pow(dxmin, 1.0/eta);
    double tmax = 1.0;
    std::cout.precision(20);    
    std::cout << "tmin set to " << tmin << std::endl;
    
//    RandomNumberGenerator g1(seed);
//    std::uniform_real_distribution<double> uniformt(tmin, 1.0);
        
    RandomNumberGenerator g2(seed+1);
    std::uniform_real_distribution<double> uniform;
    
    int ntrials = 0;
    int ngenerated = 0;
    
    std::unique_ptr<TFile> myFile( TFile::Open("NumericalInverter.root","RECREATE") );
    TH1D *h_t = new TH1D("h_t","t",100,0.0,1.0);
    TH1D *h_tp = new TH1D("h_tp","tp",100,0.027606,1.027606);
    TH1D *h_tpp = new TH1D("h_tpp","tpp",100,tmin,1.0);    
    TH1D *h_x = new TH1D("h_x","x",102,0.0,1.02);
    TH1D *h_xp = new TH1D("h_xp","xp",102,0.9900,1.0002);
    TH1D *h_xpp = new TH1D("h_xpp","xpp",102,0.99990,1.000002);
    TH1D *h_xppp = new TH1D("h_xppp","xppp",102,0.999990,1.0000002);                  

    maxweight = -1.0;
    
 // First calculation of normalization integral
    auto normintegral = Evaluate( tmin, tmax, 400, eta, a2, a3, a2p, a3p, frac );
    double normFactor = normintegral.first;
    double fracerr = (normintegral.first - normintegral.second)/normintegral.first;
    
    std::cout << "Normalization integral " << normintegral.first << " " << normintegral.second << " " << fracerr << std::endl;
    
    std::set<std::pair<double, double>> quantileSet;
    std::set<std::pair<double, double>, std::greater<std::pair<double,double>> > quantileROSet;  // reverse ordered
    
    // Construct the cumulative distribution function as function of t
    
    // First insert the end-points. Each element is the cdf value (F(t), t). 
    // So given desired quantile value of p, we want to find the value of t such that F(t) = p.
    // ie. we need to be able to invert the relationship, in the sense of measuring t = F^{-1} (p).
    quantileSet.insert(std::make_pair(0.0, tmin));
    quantileROSet.insert(std::make_pair(0.0, tmin));
    quantileSet.insert(std::make_pair(1.0, tmax));
    quantileROSet.insert(std::make_pair(1.0, tmax));        
    
    // Build-up the initial "look-up table" type container.
    // Here we sample evenly in t, and measure the values F(t) by integration.
    // The other approach would be to sample evenly in p, and compute the 
    // corresponding values of t.
    int N=1000;
    double dt = (tmax - tmin)/double(N);
    double cdf;
    for (int i = 1; i < N; i++ ) {
        double t = tmin + double(i)*dt;
        if( t < 0.285 ){
            auto integral = Evaluate( tmin, t, 400, eta, a2, a3, a2p, a3p, frac );
            cdf = integral.first/normintegral.first;
            std::cout << "i, F(t), t " << i << " " << cdf << " " << t << std::endl;
            quantileSet.insert(std::make_pair(cdf, t));
            quantileROSet.insert(std::make_pair(cdf, t));            
        }
        else if( t < 0.975 ){
            auto integral = Evaluate( t, tmax, 400, eta, a2, a3, a2p, a3p, frac );
            cdf = 1.0 - (integral.first/normintegral.first) ;           
            //std::cout << "i, F(t), t " << i << " " << cdf << " " << t << std::endl;
            quantileSet.insert(std::make_pair(cdf, t));
            quantileROSet.insert(std::make_pair(cdf, t));                       
        }    
    }
    
    // Now let's test our algorithm  
    // Choose random number (hard-coded for now) 
    
    while (ngenerated < nevents){
    
        double u = uniform(g2);
    
    std::cout << "Target quantile " << u << std::endl;
    std::pair<double, double> ptest = std::make_pair(u, 0.5); // Second element of the pair is needed - but the value should be immaterial as we have a set not a multiset.
    auto plo = quantileROSet.lower_bound(ptest);              // Get iterator to the position of the lower edge of the bracketing interval
    auto phi = quantileSet.lower_bound(ptest);                // Get iterator to the position of the upper edge of the bracketing interval 
    std::cout << " plo " << (*plo).first << " " << (*plo).second << std::endl;
    std::cout << " phi " << (*phi).first << " " << (*phi).second << std::endl;
    
    // Now we want to do some root-finding.
    // Step 1. Linear interpolation to get best estimate.
    double tlo = (*plo).second;
    double thi = (*phi).second;
    double flo = (*plo).first - u;
    double fhi = (*phi).first - u;
    
    double troot = rootFinderAlt(tlo, thi, flo, fhi, u, normFactor, tmin, tmax, eta, a2, a3, a2p, a3p, frac); 
    
    // FIXME 1 - Remember with the table saved, we only need to calculate the increments / decrements.
    // FIXME 2 - Can also save the first and second derivatives at each table abscissa.
    // FIXME 3 - Can probably calculate the second derivative efficiently and use a higher order version than NR?
    // FIXME 4 - Add a test of the integration method for a non-trivial case.
    
 // This is likely good enough. Now add these to the sets
 // Adding them to the sets may be inefficient/unwarranted after some point. We will see.
 
    int ninsertionMax = 10000;
    if(ngenerated < ninsertionMax){
        quantileSet.insert(std::make_pair(u, troot));
        quantileROSet.insert(std::make_pair(u, troot));
    }     

    double t = troot;
    ngenerated++; 
 
            h_t->Fill(t);
            h_tp->Fill(t);    // simmilar binning to ROOT histograms. 
            h_tpp->Fill(t);
   // Transform back to x
            double x = 1.0 - std::pow(t, eta);
            h_x->Fill(x);
            h_xp->Fill(x);
            h_xpp->Fill(x);
            h_xppp->Fill(x);            
            if(ngenerated <=100)std::cout << "Accept  " << ngenerated << " " << t << " " << x << " " << u <<  std::endl;
                 
    }
    
    std::cout << "Ngen " << ngenerated << std::endl;
    
    h_t->Write();
    h_tp->Write();
    h_tpp->Write();
    h_x->Write();
    h_xp->Write();
    h_xpp->Write();
    h_xppp->Write();
    myFile->Close();
 
}

int main(int argc, char **argv) {

    CLI::App app{"Model beta distribution"};
    
    int nevents = 10;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 33579L;
    app.add_option("-s,--seed", seed, "Seed for random number generator");
    
    double dxmin = 2.0e-6;
    app.add_option("-d,--dxmin", dxmin, "Cutoff parameter");
    
    double eta = 4.0;
    app.add_option("-e,--eta", eta, "eta (singularity regulation parameter)");
    
    double a2 = 14.26;
    app.add_option("--par1", a2, "a2 parameter");
    
    double a3 = -0.6233;
    app.add_option("--par2", a3, "a3 parameter");
    
    double a2p = 30.50;
    app.add_option("--par3", a2p, "a2p parameter");
    
    double a3p = -0.6721;
    app.add_option("--par4", a3p, "a3p parameter");    
    
    double frac = 0.674;
    app.add_option("-f,--frac", frac, "Fraction of first Beta component");
    
    double WTMAX = 0.52;
    app.add_option("-w,--weightmax", WTMAX, "Maximum weight");    

    CLI11_PARSE(app, argc, argv);

    std::cout << "nevents   " << nevents << std::endl;
    std::cout << "seed      " << seed << std::endl;
    std::cout << "dxmin     " << dxmin << std::endl;
    std::cout << "eta       " << eta << std::endl;
    std::cout << "a2        " << a2 << std::endl;
    std::cout << "a3        " << a3 << std::endl;  
    std::cout << "a2p       " << a2p << std::endl;  
    std::cout << "a3p       " << a3p << std::endl;                  
    std::cout << "frac      " << frac << std::endl;  
    std::cout << "WTMAX     " << WTMAX << std::endl;             
 
    UniGenerator(nevents, seed, dxmin, eta, a2, a3, a2p, a3p, frac, WTMAX);
       
    return 0;
    
}
