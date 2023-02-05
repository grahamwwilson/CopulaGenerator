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
// Generate univariate distributions of (E and mll) for A -> f1 f2 B where A and B are neutralinos
// according to equation B.76 of Baer and Tata for 3-body decays through Z couplings. 
// This equation allows for non-zero fermion mass. 
//
//                         Graham W. Wilson, 15-JAN-2023.
//

void UniGenerator(int nevents, unsigned long int seed, double mA, double mB, double mf, double signValue, double WTMAX){

// Sample from the energy distribution of neutralino B in the decay of neutralino A -> f1 f2 B
//
// Calculations are done in the rest frame of A.

    double mZ = 91.1876;
    double maxweight;
    
    double Emin = mB;
    double Emax = (mA*mA + mB*mB - 4.0*mf*mf)/(2.0*mA);
    
    RandomNumberGenerator g1(seed);
    std::uniform_real_distribution<double> uniformE(Emin, Emax);    
    RandomNumberGenerator g2(seed+1);
    std::uniform_real_distribution<double> uniform;    
        
    std::cout.precision(10);
    std::cout << "Emax  = " << Emax << std::endl;
    
    int ntrials = 0;
    int ngenerated = 0;
    
    std::unique_ptr<TFile> myFile( TFile::Open("BaerTata.root","RECREATE") );
    TH1D *h_mff = new TH1D("h_mff","Difermion Mass (GeV)",120,0.0,mA-mB);
    TH1D *h_mffsq = new TH1D("h_mffsq","Difermion Mass Squared (GeV^{2})",120,0.0,(mA-mB)*(mA-mB));
    TH1D *h_E = new TH1D("h_E","Neutralino Energy (GeV)",150,Emin,Emax);
    
    maxweight = -1.0;
    
    while (ngenerated < nevents){
    
// Try E on [Emin, Emax]
        double E = uniformE(g1);

// Evaluate unnormalized pdf
        double Bf = sqrt(1.0 - ( (4.0*mf*mf)/(mA*mA + mB*mB - 2.0*E*mB) ) );
        double densq = pow( mA*mA + mB*mB - mZ*mZ - 2.0*E*mA, 2.0);
        double prefactor = Bf*sqrt(E*E - mB*mB)/densq;
        double term1 = E*(mA*mA + mB*mB - signValue*2.0*mA*mB);
        double term2 = -mA*(E*E + mB*mB + Bf*(E*E - mB*mB)/3.0);
        double term3 = signValue*mB*(mA*mA + mB*mB - 2.0*mf*mf);
        
        double pdf = prefactor*(term1 + term2 + term3);
        
        ntrials++;
        
        if(pdf > maxweight)maxweight = pdf;
        
        double r = WTMAX*uniform(g2);       
                    
        if ( pdf > r ){
   // KEEP these random variates         
            ngenerated++; 
                
   // Transform back to masses
            double mffsq = mA*mA + mB*mB -2.0*E*mA;
            double mff = sqrt(mffsq);
            h_E->Fill(E);
            h_mffsq->Fill(mffsq);
            h_mff->Fill(mff);
            if(ngenerated <=100)std::cout << "Accept  " << ngenerated << " " 
                                << E << " " << mffsq << " " << mff<< " " << pdf <<  std::endl;
        }                 
    }
    
    std::cout << "Eigenvalue sign = " << signValue << std::endl;
    std::cout << "WTMAX set to " << WTMAX << std::endl;
    std::cout << "Max weight observed of " << maxweight << std::endl;
    if(maxweight > WTMAX)std::cout << "MAXWT exceeded. Need to increase! " << std::endl;
    std::cout << "Ntrials, Ngen, Efficiency " << ntrials << " " << ngenerated << " " << double(ngenerated)/double(ntrials) << std::endl;
    
    h_mff->Write();
    h_mffsq->Write();
    h_E->Write();
    myFile->Close();
    
}

int main(int argc, char **argv) {

    CLI::App app{"Univariate generator for 3-body Z-mediated neutralino decay, A -> f1 f2 B"};
    
    int nevents = 9;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 23579L;
    app.add_option("-s,--seed", seed, "Seed for random number generator");
    
    double mA = 300.0;
    app.add_option("-a,--ma", mA, "Parent neutralino mass (GeV)");
    
    double mB = 270.0;
    app.add_option("-b,--mb", mB, "Child neutralino mass (GeV)");    
    
    double mf = 0.0;
    app.add_option("-f,--mf", mf, "Fermion mass (GeV)");
    
    int signvalue = -1;
    app.add_option("-e,--evsign", signvalue, "Relative mass eigenvalue sign (+-1)");
    
    double WTMAX = 0.12;
    app.add_option("-w,--weightmax", WTMAX, "Maximum weight");    

    CLI11_PARSE(app, argc, argv);

    std::cout << "nevents   " << nevents << std::endl;
    std::cout << "seed      " << seed << std::endl;
    std::cout << "mA        " << mA << std::endl;
    std::cout << "mB        " << mB << std::endl;
    std::cout << "mf        " << mf << std::endl;
    double signValue = double(signvalue);
    std::cout << "signValue " << signValue << std::endl;   
    std::cout << "WTMAX     " << WTMAX << std::endl;             
 
    UniGenerator(nevents, seed, mA, mB, mf, signValue, WTMAX);
       
    return 0;
    
}
