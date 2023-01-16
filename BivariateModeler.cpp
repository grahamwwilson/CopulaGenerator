// You may choose to use #include "CLI11.hpp" to pick up the repo version header file 
#include "CLI11.hpp"
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
#include <cmath>
#include "TH2.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <map>
#include <string>

typedef std::mt19937 RandomNumberGenerator;

std::string copulafile;

//
// Generate b-variate distribution of (x,y) according to Nojiri, Yamada dG/dxdy distribution.
// Also potentially stick this through the ecdf copula stuff.
//
//                         Graham W. Wilson, 14-JAN-2023.
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

double kallen(double x, double y, double z){
// Kallen kinematical function. See Byckling and Kajantie II.6.3. 
    return std::pow(x - y - z, 2.0) - 4.0*y*z;
}

double angleFormula(double s, double mi, double mj, double sj, double sk, double si){
   //
   // Compute the opening angles in the 3-body rest frame using cyclic permutations of BK equation V.1.4.
   //                double costh12 = angleFormula(s, m1, m2, s23, s31, s12); 

       double costh;
       
       costh = ( (s + mi*mi - sj) * (s + mj*mj - sk) + 2.0*s*(mi*mi + mj*mj - si) ) /
               std::sqrt(kallen(s, mi*mi, sj) * kallen(s, mj*mj, sk) );

       return costh;
}

double decayAngleFormula(double s, double mi, double mj, double mk, double sj, double si){
   //
   // Compute the decay angle in the 2-particle rest frame using cyclic permutations of BK equation V.1.9.
   // Using the BK notation, Rij is the rest-frame for particles i and j, with pi + pj = 0.
   //
   // double costh23star = decayAngleFormula(s, m2, m3, m1, s23, s12);

       double costh;
       
       costh = ( (s - sj - mk*mk) * (sj + mi*mi - mj*mj) + 2.0*sj*(mk*mk + mi*mi - si) ) /
               std::sqrt(kallen(s, sj, mk*mk) * kallen(sj, mi*mi, mj*mj) );

       return costh;
}

void BiGenerator(int nevents, unsigned long int seed, double mA, double mB, double mf, double signValue, double WTMAX, bool passthrough){

// Sample from the 3-body phase space for A -> l1 l2 B based on the presentation 
// in Nojiri and Yamada. This currentyl neglects fermion masses.
// Althought at least it should be relatively trivial to incorporate 
// fermion masses on the allowed phase-space if not the 2-d matrix 
// element-squared itself.
//
// x = [m(l1,N)/mA]**2
// y = [m(l2,N)/mA]**2
// z = [m(l1,l2)/mA]**2
// 
// We define r = +- (mB/mA) where both signs are feasible
// rz = mZ/mA
// Could it be that rZ can also be negative?? No it doesn't matter.

    double mZ = 91.1876;
    double pdfave = 1.25391;     //for signValue = -1
//    double pdfave=3.23978;     // for signValue = +1
    double maxweight = -1.0;
    double rB = signValue*(mB/mA);
    double rBsq = rB*rB;
    double rZsq = mZ*mZ/(mA*mA);
    double zsqmax = std::pow( (mA-mB)/mA, 2.0);

    RandomNumberGenerator g1(seed);
    std::uniform_real_distribution<double> uniformr(rBsq, 1.0);    
    RandomNumberGenerator g2(seed+1);
    std::uniform_real_distribution<double> uniform;    
    
    std::set<std::pair<std::pair<double,double>, int>> uvset;
    std::set<std::pair<std::pair<double,double>, int>> vuset; 
    
    double u1,u2;
    
    int ntrials = 0;
    int ngenerated = 0;
    int noutofbounds = 0;
    int ninbounds = 0;
    
    double testQuantity;
    
    std::unique_ptr<TFile> myFile( TFile::Open("PhaseSpace.root","RECREATE") );
    TH1D *h_mll = new TH1D("h_mll","Dilepton Mass (GeV)",120,0.0,mA-mB);
    TH1D *h_mllR = new TH1D("h_mllR","Reweighted Dilepton Mass (GeV)",120,0.0,mA-mB);    
    TH1D *h_costh12 = new TH1D("h_costh12","costh12",100,-1.0,1.0);
    TH1D *h_costh23 = new TH1D("h_costh23","costh23",100,-1.0,1.0);
    TH1D *h_costh31 = new TH1D("h_costh31","costh31",100,-1.0,1.0); 
    TH1D *h_pt = new TH1D("h_pt","pT (lepton) (GeV)",100,0.0,0.5*(mA-mB));
    TH1D *h_costh12star = new TH1D("h_costh12star","costh12star",100,-1.0,1.0);
    TH1D *h_costh12starR = new TH1D("h_costh12starR","Reweighted costh12star",100,-1.0,1.0);    
    TH1D *h_costh23star = new TH1D("h_costh23star","costh23star",100,-1.0,1.0); 
    TH1D *h_costh31star = new TH1D("h_costh31star","costh31star",100,-1.0,1.0); 
    TH2D *h_xy = new TH2D("h_xy","y vs x",100,rBsq,1.0,100,rBsq,1.0); 
    TH2D *h_xz = new TH2D("h_xz","x vs z",100,0.0,zsqmax,100,rBsq,1.0);
    TH2D *h_xzR = new TH2D("h_xzR","Reweighted x vs z",100,0.0,zsqmax,100,rBsq,1.0);        
    
    double pdfsum = 0.0;                           
    
    while (ngenerated < nevents){
    
// Try (x,y) on [rB**2, 1.0] x [rB**2, 1.0]
        double x = uniformr(g1);
        double y = uniformr(g1);
        double z = 1.0 + rBsq - x - y;
        testQuantity = z*(x*y - rBsq);
        
        ntrials++;
        
        if (ntrials <= 100){
           std::cout << "trial " << ntrials << " " << x << " " << y << " " << z << " " << testQuantity << std::endl;
        }        
        
        if (testQuantity >= 0.0) {
  // We're in the allowed region of the Dalitz plot. Need to evaluate unnormalized pdf and throw a proper random number
            ninbounds++;
            double pdf = ( (1.0-x)*(x-rBsq) + (1.0-y)*(y-rBsq) + 2.0*rB*z ) / std::pow(z - rZsq, 2.0);
            if(pdf > maxweight)maxweight = pdf;
            
            // std::cout << pdf << " x, y, z " << x << " " << y << " " << z << " pdf = " << pdf << std::endl;
            
            double r = WTMAX*uniform(g2);
            
            if ( passthrough || pdf > r ){
   // KEEP these random variates         
                ngenerated++; 
                
                pdfsum += pdf;
                 
   // Transform back to masses
                double mll  = sqrt(z*mA*mA);
                double ml1N = sqrt(x*mA*mA);
                double ml2N = sqrt(y*mA*mA);
                if(ngenerated <= 100)std::cout << "Accept  " << mll << " " << ml1N << " " << ml2N << std::endl; 
                
   // Set up variables as in Byckling and Kajantie chapter 5. But let's dispense with the confusing s1,s2,s3
                double s   = mA*mA;
                double s12 = mll*mll;          // Was s1
                double s23 = ml2N*ml2N;        // Was s2
                double s31 = ml1N*ml1N;        // Was s3
                double m1  = mf;
                double m2  = mf;
                double m3  = mB;
                
   // Compute the opening angles in the 3-body rest frame using cyclic permutations of equation V.1.4.
                double costh12 = angleFormula(s, m1, m2, s23, s31, s12);
                double costh23 = angleFormula(s, m2, m3, s31, s12, s23);
                double costh31 = angleFormula(s, m3, m1, s12, s23, s31);
   // Compute the energies and momenta in the 3-body rest frame             
                double E1 = (s + m1*m1 - s23)/(2.0*std::sqrt(s));
                double E2 = (s + m2*m2 - s31)/(2.0*std::sqrt(s));
                double E3 = (s + m3*m3 - s12)/(2.0*std::sqrt(s));             
                double p1 = std::sqrt(E1*E1 - m1*m1);
                double p2 = std::sqrt(E2*E2 - m2*m2);
                double p3 = std::sqrt(E3*E3 - m3*m3);
   // Compute the "scattering angle" in the two-body rest frame using equation V.1.9             
                double costh23star = decayAngleFormula(s, m2, m3, m1, s23, s12);   //costh_12^R23
                double costh31star = decayAngleFormula(s, m3, m1, m2, s31, s23);   //costh_23^R31
                double costh12star = decayAngleFormula(s, m1, m2, m3, s12, s31);   //costh_31^R12
                
   // Lepton pt with respect to neutralino3 direction
                double pt = p1*sin(acos(costh31));                                 

//        fout << std::setw(10) << i+1 << "  " << std::scientific << u1 << "  " << std::scientific << u2 << "  " << idistbn << std::endl;

 // At this stage (u1,u2) are random variables on (0,1) from the copula distribution

                std::pair<double,double> puv=std::make_pair(x,y);
                uvset.insert(std::make_pair(puv,ngenerated));              // Will be sorted by increasing values of u1
                std::pair<double,double> pvu=std::make_pair(y,x);
                vuset.insert(std::make_pair(pvu,ngenerated));              // Will be sorted by increasing values of u2
                
                h_mll->Fill(mll);
                h_costh12->Fill(costh12);
                h_costh23->Fill(costh23);
                h_costh31->Fill(costh31);
                h_pt->Fill(pt);
                h_costh12star->Fill(costh12star); 
                h_costh23star->Fill(costh23star); 
                h_costh31star->Fill(costh31star);
                h_xy->Fill(x,y);
                h_xz->Fill(z,x);  
                
                // reweighting tests under either Odd or Even hypothesis
                double wt = pdf/pdfave;
                
                h_mllR->Fill(mll,wt);
                h_costh12starR->Fill(costh12star,wt);
                h_xzR->Fill(z,x,wt);
                
            }
        }    
        else{ 
            noutofbounds++;
        }
    }
    
    std::cout << "rBsq " << rBsq << " sign = " << signValue << std::endl;
    std::cout << "WTMAX set to " << WTMAX << std::endl;
    std::cout << "Max weight observed of " << maxweight << std::endl;
    std::cout << "Ntrials, NOB, NIB, Ngen, Efficiency, Efficiency2 " << ntrials << " " << noutofbounds << " " << ninbounds 
              << " " << ngenerated << " " << double(ngenerated)/double(ntrials) << " " << double(ngenerated)/double(ninbounds) << std::endl;
    std::cout << "Mean pdf value for generated events " << pdfsum/double(ngenerated) << std::endl;
    
    h_mll->Write();
    h_costh12->Write();
    h_costh23->Write();
    h_costh31->Write();
    h_pt->Write();
    h_costh12star->Write();
    h_costh23star->Write();
    h_costh31star->Write();
    h_xy->Write();
    h_xz->Write();
    h_mllR->Write();
    h_costh12starR->Write();
    h_xzR->Write();
    myFile->Close();
    
// Process the sets and write out the empirical copula file   
//    Formatter(uvset, vuset);
    
}

int main(int argc, char **argv) {

    CLI::App app{"Sample from 3-body decay phase-space with appropriate matrix element"};
    
    int nevents = 1000000;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 13579L;
    app.add_option("-s,--seed", seed, "Seed");
    
    double mA = 300.0;
    app.add_option("-a,--ma", mA, "Parent neutralino mass (GeV)");
    
    double mB = 270.0;
    app.add_option("-b,--mb", mB, "Child neutralino mass (GeV)");    
    
    double mf = 0.0;
    app.add_option("-f,--mf", mf, "Fermion mass (GeV)");
    
    int signvalue = -1;
    app.add_option("-e,--evsign", signvalue, "Relative mass eigenvalue sign (+-1)");
    
    double WTMAX = 6.0;
    app.add_option("-w,--weightmax", WTMAX, "Maximum weight");     

    std::string filename = "CopulaGen.EDAT";
    app.add_option("-o,--outputfile", filename, "Output copula file");
    
    bool passthrough = false;
    app.add_flag("-p,--passthrough", passthrough, "Disable matrix element weighting");    
    
    CLI11_PARSE(app, argc, argv);

    std::cout << "nevents     " << nevents << std::endl;
    std::cout << "seed        " << seed << std::endl;
    std::cout << "mA          " << mA << std::endl;
    std::cout << "mB          " << mB << std::endl;
    std::cout << "mf          " << mf << std::endl;
    double signValue = double(signvalue);
    std::cout << "signValue   " << signValue << std::endl;   
    std::cout << "WTMAX       " << WTMAX << std::endl;    
    std::cout << "filename    " << filename << std::endl;
    std::cout << "passthrough " << passthrough << std::endl;    

    copulafile = filename;             // Write filename to globally declared copulafile string
    
    BiGenerator(nevents, seed, mA, mB, mf, signValue, WTMAX, passthrough);    
       
    return 0;
    
}
