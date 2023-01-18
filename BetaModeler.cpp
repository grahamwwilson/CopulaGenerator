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

void UniGenerator(int nevents, unsigned long int seed, double dxmin, double eta, double a2, double a3, double a2p, double a3p, double frac, double WTMAX){

    double maxweight;
    
    double tmin = std::pow(dxmin, 1.0/eta);
    std::cout.precision(10);    
    std::cout << "tmin set to " << tmin << std::endl;
    
    RandomNumberGenerator g1(seed);
    std::uniform_real_distribution<double> uniformt(tmin, 1.0);    
    RandomNumberGenerator g2(seed+1);
    std::uniform_real_distribution<double> uniform;
    
    int ntrials = 0;
    int ngenerated = 0;
    
    std::unique_ptr<TFile> myFile( TFile::Open("BetaTest.root","RECREATE") );
    TH1D *h_t = new TH1D("h_t","t",100,0.0,1.0);
    TH1D *h_tp = new TH1D("h_tp","tp",100,0.027606,1.027606);
    TH1D *h_tpp = new TH1D("h_tpp","tpp",100,tmin,1.0);    
    TH1D *h_x = new TH1D("h_x","x",102,0.0,1.02);
    TH1D *h_xp = new TH1D("h_xp","xp",102,0.9900,1.0002);
    TH1D *h_xpp = new TH1D("h_xpp","xpp",102,0.99990,1.000002);
    TH1D *h_xppp = new TH1D("h_xppp","xppp",102,0.999990,1.0000002);                  

    maxweight = -1.0;
    
    while (ngenerated < nevents){
    
// Try E on [Emin, Emax]
        double t = uniformt(g1);

// Evaluate unnormalized pdf
 
        double beta1 = std::pow(1.0 - std::pow(t,eta), a2) * std::pow(t, eta*(1.0+a3) - 1.0);
        double beta2 = std::pow(1.0 - std::pow(t,eta), a2p) * std::pow(t, eta*(1.0+a3p) - 1.0);
        double pdf = frac*beta1 + (1.0-frac)*beta2;        
        
        ntrials++;
        
        if(pdf > maxweight)maxweight = pdf;
        
        double r = WTMAX*uniform(g2);       
                    
        if ( pdf > r ){
   // KEEP these random variates         
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
            if(ngenerated <=100)std::cout << "Accept  " << ngenerated << " " << t << " " << x << " " << pdf <<  std::endl;
        }                 
    }
    
    std::cout << "WTMAX set to " << WTMAX << std::endl;
    std::cout << "Max weight observed of " << maxweight << std::endl;
    if(maxweight > WTMAX)std::cout << "MAXWT exceeded. Need to increase! " << std::endl;
    std::cout << "Ntrials, Ngen, Efficiency " << ntrials << " " << ngenerated << " " << double(ngenerated)/double(ntrials) << std::endl;
    
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
    
    int nevents = 1000;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 23579L;
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
