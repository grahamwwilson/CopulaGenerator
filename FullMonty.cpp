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
// Generate all four parts of the luminosity spectrum.
//
//                         Graham W. Wilson, 17-JAN-2023.
//

bool generateArmDoubleBeta(double t, double r){
#include "armsDoubleBeta.h"
// Evaluate unnormalized pdf
    double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
    double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
    double pdf = frac*beta1 + (1.0-frac)*beta2;
    
    bool success = false;
    
    if ( pdf > WTMAX ){
        std::cout << "Warning maximum allowed weight exceeded. Need to adjust " << pdf << " " << WTMAX << std::endl;
    }
        
    if ( pdf > WTMAX*r ) success = true;
    
    return success;
   
}

double functionGradient(double t, double normFactor){
// 
// Evaluate the function gradient at t
//
#include "bodyDoubleBeta.h" 

    double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
    double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
    double tintegrand = frac*beta1 + (1.0-frac)*beta2;
    double pdf = tintegrand/normFactor;          

    return pdf;
    
}

std::pair<double, double> Evaluate(double tmin, double tmax, int Nbase){
// 
// Evaluate the un-normalized cumulative distribution function using Boole's rule.
// This is a 5-point closed Newton-Cotes formula. 
// (the first call will be to calculate the integral over the full range)
//
#include "bodyDoubleBeta.h"

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

double rootFinderAlt( double tlo, double thi, double flo, double fhi, double u, double normFactor, double tmin, double tmax){
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
            dintegral = Evaluate( troot, thi, Nsteps );
            cdfestimate = (fhi + u) - (dintegral.first/normFactor);
        }
        else{
            dintegral = Evaluate( tlo, troot, Nsteps );
            cdfestimate = (flo + u) + (dintegral.first/normFactor);            
        }    
            
        double fn = cdfestimate - u;
        double fngrad = functionGradient( troot, normFactor );
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

std::pair<double, double> copulaVariates(double which, double u1, double v2){

// Hard code the best fit parameters based on fit to GP4X5X data (583,584 events).
// Fit parameters 
// Clayton dependence parameter 0.04969906094
// AMH dependence parameter     0.3598358846
// Clayton phiwt parameter      0.4436183234 wt = 0.8157784722

    double thetaC  = 0.04969906094;
    double thA     = 0.3598358846;
    double weightC = 0.8157784722; 
    
    double u2;
    
    if ( which <= weightC ){
  // Clayton (also known by other names. Joe (2014) calls it the bivariate Mardia-Takahasi-Clayton-Cook-Johnson copula
         
         u2 = pow(1.0 - pow(u1, -thetaC) + pow( 1.0/(v2*pow(u1, 1.0+thetaC)), thetaC/(1.0+thetaC) ) , -1.0/thetaC);      
    }
    else{
  // Ali-Mikhail-Haq
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
    // u1, u2 are the copula parameters, but since we use these to look up values in the t distribution we should "survival" them
    return std::make_pair(1.0-u1, 1.0-u2);
}


void FullMontyGenerator(int nevents, unsigned long int seed, double dxmin, double eta, double dmu1, double dmu2, double sigma1, double sigma2){

#include "regionProbabilities.h"

    double maxweight;
    
    double tmin = std::pow(dxmin, 1.0/eta);
    double tmax = 1.0;
    std::cout.precision(20);    
    std::cout << "tmin set to " << tmin << std::endl;
    
// Currently designed such that all the random number generation is done in this function itself 
// and managed here.
// Only the random numbers rather than the generators need to be passed to the called functions.

    RandomNumberGenerator gent(seed);
    std::uniform_real_distribution<double> uniformt(tmin, 1.0);
        
    RandomNumberGenerator genu(seed+1);
    std::uniform_real_distribution<double> uniform;
    
    RandomNumberGenerator geng(seed+2);
    std::normal_distribution<double> gaussian(0.0, 1.0);  // Make sure this is standardized
    
    int ntrials = 0;
    int ngenerated = 0;
    
    std::unique_ptr<TFile> myFile( TFile::Open("FullMonty.root","RECREATE") );
    TH1D *h_t = new TH1D("h_t","t",100,0.0,1.0);
    TH1D *h_tp = new TH1D("h_tp","tp",100,0.027606,1.027606);
    TH1D *h_tpp = new TH1D("h_tpp","tpp",100,tmin,1.0);    
    TH1D *h_x = new TH1D("h_x","x",102,0.0,1.02);
    TH1D *h_xp = new TH1D("h_xp","xp",102,0.9900,1.0002);
    TH1D *h_xpp = new TH1D("h_xpp","xpp",102,0.99990,1.000002);
    TH1D *h_xppp = new TH1D("h_xppp","xppp",102,0.999990,1.0000002);
    TH1D *h_ecm = new TH1D("h_ecm","ecm",11000,0.0,1.10);                     

    maxweight = -1.0;
    
 // First calculation of normalization integral
    auto normintegral = Evaluate( tmin, tmax, 400);
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
            auto integral = Evaluate( tmin, t, 400);
            cdf = integral.first/normintegral.first;
            std::cout << "i, F(t), t " << i << " " << cdf << " " << t << std::endl;
            quantileSet.insert(std::make_pair(cdf, t));
            quantileROSet.insert(std::make_pair(cdf, t));            
        }
        else if( t < 0.975 ){
            auto integral = Evaluate( t, tmax, 400);
            cdf = 1.0 - (integral.first/normintegral.first) ;           
            //std::cout << "i, F(t), t " << i << " " << cdf << " " << t << std::endl;
            quantileSet.insert(std::make_pair(cdf, t));
            quantileROSet.insert(std::make_pair(cdf, t));                       
        }    
    }
    
    // Now let's test our algorithm  
    // Choose random number (hard-coded for now) 
    
    double x1,x2;     // These are the unsmeared quantities we want!
    double x1p,x2p;   // These are the beam energy spread smeared quantities.
    
    std::string region;
    
    while (ngenerated < nevents){
    
        double uregion = uniform(genu);
        
        if( uregion < ppeak ){
            // Peak region
            std::cout << "PEAK " << std::endl;
            region = "PEAK "; 
            x1 = 1.0;
            x2 = 1.0;
        }
        else if ( uregion < ppeak + parm1 ){
            // Arm1
            // Generate arm specific double beta distribution for x1
            std::cout << "ARM1 " << std::endl;
            region = "ARM1 ";             
            bool arm1Failure = true;
            double t1, rtest;
            while (arm1Failure) {
                t1 = uniformt(gent);
                rtest = uniform(genu);
                bool success = generateArmDoubleBeta(t1, rtest);
                if(success){
                   x1 = 1.0 - std::pow(t1, eta);
                   arm1Failure = false;
                }
            }
            x2 = 1.0;
        }
        else if ( uregion < ppeak + parm1 + parm2 ){
            // Arm2 
            // Generate same arm specific double beta distribution for x2
            std::cout << "ARM2 " << std::endl;
            region = "ARM2";                
            bool arm2Failure = true;
            double t2, rtest;
            while (arm2Failure) {
                t2 = uniformt(gent);
                rtest = uniform(genu);
                bool success = generateArmDoubleBeta(t2, rtest);
                if(success){
                   x2 = 1.0 - std::pow(t2, eta);
                   arm2Failure = false;
                }
            }            
            x1 = 1.0;
        }
        else{
            // Body
            std::cout << "BODY " << std::endl; 
            region = "BODY ";            
            // First generate copula (u,v) values
            double which = uniform(genu);
            double u1 = uniform(genu);
            double v2 = uniform(genu);
            std::pair<double,double> cpair = copulaVariates(which, u1, v2);
            double u = cpair.first;
            double v = cpair.second;
    
    // First variable
            std::cout << "Target quantile " << u << std::endl;
            std::pair<double, double> ptest = std::make_pair(u, 0.5); // Second element of the pair is needed - but the value should be immaterial as we have a set not a multiset.
            auto plo = quantileROSet.lower_bound(ptest);              // Get iterator to the position of the lower edge of the bracketing interval
            auto phi = quantileSet.lower_bound(ptest);                // Get iterator to the position of the upper edge of the bracketing interval 
            std::cout << " plo " << (*plo).first << " " << (*plo).second << std::endl;
            std::cout << " phi " << (*phi).first << " " << (*phi).second << std::endl;
    
    // Now we want to do some root-finding.
            double tlo = (*plo).second;
            double thi = (*phi).second;
            double flo = (*plo).first - u;
            double fhi = (*phi).first - u;
            double troot = rootFinderAlt(tlo, thi, flo, fhi, u, normFactor, tmin, tmax); 
    
            int ninsertionMax = 10000;
            if(ngenerated < ninsertionMax){
                quantileSet.insert(std::make_pair(u, troot));
                quantileROSet.insert(std::make_pair(u, troot));
            }     

            double t1 = troot;
            x1 = 1.0 - std::pow(t1, eta);
            
    // Second variable        
            std::cout << "Target quantile " << v << std::endl;
            ptest = std::make_pair(v, 0.5); // Second element of the pair is needed - but the value should be immaterial as we have a set not a multiset.
            plo = quantileROSet.lower_bound(ptest);              // Get iterator to the position of the lower edge of the bracketing interval
            phi = quantileSet.lower_bound(ptest);                // Get iterator to the position of the upper edge of the bracketing interval 
            std::cout << " plo " << (*plo).first << " " << (*plo).second << std::endl;
            std::cout << " phi " << (*phi).first << " " << (*phi).second << std::endl;
    
    // Now we want to do some root-finding.
            tlo = (*plo).second;
            thi = (*phi).second;
            flo = (*plo).first - v;
            fhi = (*phi).first - v;
            troot = rootFinderAlt(tlo, thi, flo, fhi, v, normFactor, tmin, tmax); 
    
            if(ngenerated < ninsertionMax){
                quantileSet.insert(std::make_pair(v, troot));
                quantileROSet.insert(std::make_pair(v, troot));
            }     

            double t2 = troot;
            x2 = 1.0 - std::pow(t2, eta);            
            
       }
       ngenerated++;
       
       // Gaussian smearing associated with initial beam
         
       x1p = x1 + sigma1*gaussian(geng);
       x2p = x2 + sigma2*gaussian(geng);
       double sqrts = std::sqrt(x1p*x2p);     
       h_ecm->Fill(sqrts);
            
       if(ngenerated <=100)std::cout << "Accept  " << ngenerated << " " << region << " (x1, x2) " << x1 << " " << x2 <<  std::endl;
       
/* 
            h_t->Fill(t);
            h_tp->Fill(t);    // simmilar binning to ROOT histograms. 
            h_tpp->Fill(t);
   // Transform back to x
            double x = 1.0 - std::pow(t, eta);
            h_x->Fill(x);
            h_xp->Fill(x);
            h_xpp->Fill(x);
            h_xppp->Fill(x);  
*/                   
                 
    }
    
    std::cout << "Ngen " << ngenerated << std::endl;
    
    h_t->Write();
    h_tp->Write();
    h_tpp->Write();
    h_x->Write();
    h_xp->Write();
    h_xpp->Write();
    h_xppp->Write();
    h_ecm->Write();
    myFile->Close();
 
}

int main(int argc, char **argv) {

    CLI::App app{"Parametrized simulation of full luminosity spectrum"};
    
    int nevents = 10;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 33579L;
    app.add_option("-s,--seed", seed, "Seed for random number generator");
    
    double dxmin = 2.0e-6;
    app.add_option("-d,--dxmin", dxmin, "Cutoff parameter");
    
    double eta = 4.0;
    app.add_option("-e,--eta", eta, "eta (singularity regulation parameter)");
       
// Now hard-coded in armsDoubleBeta.h where it is used.
//    double WTMAX = 0.52;
//    app.add_option("-w,--weightmax", WTMAX, "Maximum weight");

    double dmu1 = 0.0;
    app.add_option("--dmu1", dmu1, "Electron energy scale deviation (ppm)"); 
    
    double dmu2 = 0.0;
    app.add_option("--dmu2", dmu2, "Positron energy scale deviation (ppm)");  

    double sigma1 = 0.190e-2;
    app.add_option("--sigma1", sigma1, "Electron beam energy spread"); 
    
    double sigma2 = 0.152e-2;
    app.add_option("--sigma2", sigma2, "Positron beam energy spread");        
    

    CLI11_PARSE(app, argc, argv);

    std::cout << "nevents   " << nevents << std::endl;
    std::cout << "seed      " << seed << std::endl;
    std::cout << "dxmin     " << dxmin << std::endl;
    std::cout << "eta       " << eta << std::endl;
    std::cout << "dmu1      " << dmu1 << std::endl;
    std::cout << "dmu2      " << dmu2 << std::endl;     
    std::cout << "sigma1    " << sigma1 << std::endl;
    std::cout << "sigma2    " << sigma2 << std::endl;             
 
    FullMontyGenerator(nevents, seed, dxmin, eta, dmu1, dmu2, sigma1, sigma2);
       
    return 0;
    
}
