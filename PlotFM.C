void PlotFM(string histfile="h2", string plot="Plot.png", double ymin = 2000.0, double ymax = 80000.0){

    gStyle->SetOptStat(" ");
    gStyle->SetLegendTextSize(0.04);
    
    TCanvas *c1 = new TCanvas("c1","multipads",1600,1200);
    
    TFile *f1 = new TFile("FullMonty-2000000-C0.root");
    TFile *f2 = new TFile("FullMonty-2000000-C2.root");
    TFile *f3 = new TFile("FullMonty-2000000-C1.root");        

    c1->SetTicks(1,1);
    c1->SetGrid(1,1);
    c1->SetLeftMargin(0.15);
  
    TH1D *hSS = (TH1D*) f1->Get(histfile.c_str());
    TH1D *hOS = (TH1D*) f2->Get(histfile.c_str());
    TH1D *hPS = (TH1D*) f3->Get(histfile.c_str());      
    
    hSS->SetMaximum(ymax);
    hSS->SetMinimum(ymin);
    
    hSS->SetTitle("ILC 250-SetA Beam Parameters");       
    
    hSS->GetXaxis()->SetTitle("Center-of-Mass Energy (GeV)");
    hSS->GetYaxis()->SetTitle("Events per 50 MeV bin");
    hSS->SetLineWidth(3);   
    hSS->SetLineColor(kBlue);

    hSS->Draw("e");
    hOS->SetLineColor(kRed);
    hOS->SetLineWidth(3);
    hOS->Draw("histsame");
    
    hPS->SetLineColor(kViolet);
    hPS->SetLineWidth(3);
    //hPS->Draw("histsame");
    
    hSS->Draw("esame");
        
    TLegend* leg = new TLegend(0.62, 0.66, 0.62 + 0.30, 0.66 + 0.20);
    leg->SetTextFont(42);
       
//    leg->SetHeader("ME models","C");                 

    TLegendEntry *entry=leg->AddEntry("hSS","No copula","fl");
    entry->SetLineColor(kBlue);
    entry->SetLineWidth(3);
              
    TLegendEntry *entry3=leg->AddEntry("hOS","Fit copula","fl");
    entry3->SetLineColor(kRed);
    entry3->SetLineWidth(3);
    
    TLegendEntry *entry2=leg->AddEntry("hPS","Independence copula","fl");
    entry2->SetLineColor(kViolet);
    entry2->SetLineWidth(3);    
                
    leg->Draw();
 
    TText *t = new TText(0.5,-0.0025,"Clayton + Ali-Mikhail-Haq Mixture Model (3 parameters)");
//    TText *t = new TText(0.5,-0.0025,"Ali-Mikhail-Haq Model (1 parameter)"); 
//    TText *t = new TText(0.5,-0.0025,"Clayton Model (1 parameter)");       
    t->SetTextAlign(22);
    t->SetTextColor(kRed+2);
    t->SetTextFont(43);
    t->SetTextSize(40);
//    t->SetNDC();
//    t->Draw();       
    
    c1->Print("Plot.png");

}
