void example() {
  //Open root file 
  TFile *DATA = TFile::Open("/home/snuintern1/skflat/SKFlatAnalyzer/SSdimuon_Analysis/DATA/SSdimuon_SingleMuon_2016.root");

  TCanvas *c1 = new TCanvas("c1","c1");
  TH1F* h1 = new TH1F();
  TH1F* h2 = new TH1F(); 
 
  TString a,b,c;

  vector<TString> IDs;
  IDs = {"POGLoose_Central","POGTight_Central"};

  vector<TString> Iso;
  Iso = {"","_Iso","_Iso_NonIso"};

  vector<TString> Jet;
  Jet = {"","_1j","_2j","_MET40","_MET40_1j","_MET40_2j","_MET40_1b","MET40_1j1b","MET40_2j1b"};

        a = "DATA/"+IDs.at(0)+"/Mll_pp"+Iso.at(0)+Jet.at(0); cout << a << endl;
        b = "DATA/"+IDs.at(0)+"/Mll_mm"+Iso.at(0)+Jet.at(0); cout << b << endl;
        c = "DATA_"+IDs.at(0)+"_Ratio"+Iso.at(0)+Jet.at(0)+".pdf"; cout << c << endl;

        h1 = (TH1F*)DATA->Get(a);
        h2 = (TH1F*)DATA->Get(b);

        h1->Draw("E P");
        h1->Rebin(10);
        h1->GetXaxis()->SetRangeUser(10.,3000.);
        h1->SetMarkerSize(8);
        h1->SetMarkerColor(40);
        h1->Sumw2();
        h1->GetXaxis()->SetTitle("M(ll) [GeV]");
        h1->GetYaxis()->SetTitle("Events/10GeV");
        
        h2->Draw("E P SAME"); 
        h2->Rebin(10);
        h2->GetXaxis()->SetRangeUser(10.,3000.);
        h2->SetMarkerSize(8);
        h2->SetMarkerColor(30);
        h2->Sumw2();
       
        auto rp = new TRatioPlot(h1,h2,"divsym");
        rp->Draw();
        rp->GetUpperPad()->SetLogy();
        rp->GetLowerRefYaxis()->SetTitle("Ratio(++/--)");
        rp->GetLowerRefYaxis()->SetRangeUser(0.5,1.5);
        rp->SetH1DrawOpt("E P");

        c1->Update();
        
        c1->SaveAs(c);
        h1=NULL;
        h2=NULL;
        c1->Clear();
}
