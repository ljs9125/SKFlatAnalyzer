#include <iostream>
#include <string>

void SSlepton_3() {

  gStyle->SetOptLogy(1); 
  gStyle->SetOptLogx(1); 

  //Prepare for making directory name
  std::string list[] ={"DATA","WZ_pythia","ZZ_pythia","DYJets"};

  std::string a, b, c;
  a =  "/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/";
  b =  ".root";

  //Prepare for importing histogram
  std::string type[] = {"SS", "OS", "++", "--"};
  std::string d, e;
  d =  "POGMedium_Central";
  
  //Make Canvas 
  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  TCanvas *c3 = new TCanvas("c3","c3");
  TCanvas *c4 = new TCanvas("c4","c4");
 
  THStack *hs = new THStack("hs","");

  c2->Divide(2,2);

  for (int i=0; i<4;i++) {

    //Make directory name + histogram
    if (list[i]=="DATA") {
      c = a + "DATA/SSlepton_SingleMuon.root";
      cout << c << endl;
    }
    else {
      c = a + "SSlepton_"+ list[i] + b;
      cout << c << endl;
    }

    char f[c.size() + 1]; 
    c.copy(f, c.size()+1);
    f[c.size()]='\0';

   
    //Open root file from directory
    TFile *file = TFile::Open(f);

    for (int j=0; j < 4; j++) {
      e  = d + "/" + type[j] + "_" + list[i] + "_" + d;
      cout << e << endl;

      char h[e.size()+1];
      e.copy(h,e.size()+1);
      h[e.size()]='\0';

      TH1F *hi = (TH1F*)file->Get(h);

      if (type[j]=="SS"){
         c1->cd();
         hi->Rebin(10);
         hi->GetXaxis()->SetRangeUser(10.,3000.);
         if (list[i]=="DATA"){
           hi->Draw("E P");
           hi->SetMarkerStyle(8);
         }
         else if (list[i] == "WZ_pythia"|list[i] == "ZZ_pythia"){
           hi->Draw("HIST SAME");
           hi->SetFillColor(kRed);
           hs->Add(hi);
           
         }
         else {
           hi->Draw("HIST SAME");
           hi->SetFillColor(kGreen);
           hs->Add(hi);
         }
         hs->Draw("HIST SAME");
      }
      else {
        c2->cd(j+1);
        if (list[i]=="DATA"){ 
          hi->Rebin(10);
          hi->GetXaxis()->SetRangeUser(10.,3000.);
          hi->Draw("E P");
          hi->SetMarkerStyle(8);
        }
      }

    } 

  }
/*
  TH1F *h1 = (TH1F*)DATA->Get("POGMedium_Central/SS_DATA_POGMedium_Central");
  TH1F *h2 = (TH1F*)WZ_pythia->Get("POGMedium_Central/SS_WZ_pythia_POGMedium_Central");
  TH1F *h3 = (TH1F*)ZZ_pythia->Get("POGMedium_Central/SS_ZZ_pythia_POGMedium_Central");
  TH1F *h4 = (TH1F*)DYJets->Get("POGMedium_Central/SS_DYJets_POGMedium_Central");
 

// THStack *hs = new THStack("hs",""); 


  //draw data 
  c->cd();

  h1->Rebin(10);
  h1->Draw("E P");
  h1->GetXaxis()->SetRangeUser(10.,3000.);
  h1->SetMarkerStyle(8);

  //draw MC Samples 1)charge mismatch 2)misidentification 3)prompt SS
  //Using THStack for MC Sample to compare with data
  
  
  h3->Add(h2);

  h3->Rebin(10);
  h3->Draw("HIST SAME");
  h3->SetLineColor(kGreen);
  h3->SetFillStyle(3003);
  h3->SetFillColor(kGreen);
  h3->GetXaxis()->SetRangeUser(10.,3000.);

  hs->Add(h3);
   
  h4->Rebin(10);
  //h4->Draw("HIST SAME");
  h4->SetLineColor(kRed);
  h4->SetFillStyle(3003);
  h4->SetFillColor(kRed);

  hs->Add(h4);

  hs->Draw("HIST SAME");

  
  //Axis Label
  hi->GetYaxis()->SetTitle("Entries");
  hi->GetXaxis()->SetTitle("M(mu+mu-) [GeV]");
  hi->SetTitle("SameSign dimuon invariant mass distribution"); 

  //Legend
  TLegend *l = new TLegend(0.2,0.3,0.4,0.4);
  l->AddEntry(hi,"DATA","lep"); 
  l->Draw();
  
  //Save as pdf file
  c->SaveAs("./output/SSdimuon.pdf");
*/
}
