void SSdimuon_150to300() {
  //Open root file 
  TFile* DATA = TFile::Open("/home/snuintern1/skflat/SKFlatAnalyzer/SSdimuon_Analysis/DATA/SSdimuon_SingleMuon_2016.root");

  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F* h1 = new TH1F();
  TH1F* h2 = new TH1F();  
 
  TString a,b,c;

  vector<TString> IDs;
  IDs = {"POGLoose_Central","POGTight_Central"};

  vector<TString> Iso;
  Iso = {"","_Iso","_Iso_NonIso"};

  vector<TString> Jet;
  Jet = {"","_1j","_2j","_MET40","_MET40_1j","_MET40_2j","_MET40_1b","_MET40_1j1b","_MET40_2j1b"};

  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<9;k++){
        a = "DATA/"+IDs.at(i)+"/Mll_pp"+Iso.at(j)+Jet.at(k);
        b = "DATA/"+IDs.at(i)+"/Mll_mm"+Iso.at(j)+Jet.at(k);
        c = "hist/DATA/150to300/"+IDs.at(i)+"/Ratio"+Iso.at(j)+Jet.at(k)+".pdf";

        h1 = (TH1F*)DATA->Get(a);
        h2 = (TH1F*)DATA->Get(b);

        c1->cd();

        h1->Draw("E P");
        h1->GetXaxis()->SetRangeUser(150.,300.);
        h1->GetXaxis()->SetTitle("M(ll) [GeV]");
        h1->GetYaxis()->SetTitle("Events/10GeV");
        h1->SetMarkerSize(7);
        h1->Sumw2();

        h2->Draw("SAME E P");
        h2->GetXaxis()->SetRangeUser(150.,300.);
        h2->SetMarkerSize(7);
        h2->SetMarkerColor(2);
        h2->SetLineColor(2);

        auto rp = new TRatioPlot(h1,h2,"divsym");
        rp->SetH1DrawOpt("E P");
        rp->SetH2DrawOpt("E P");
        rp->Draw("confint");
        rp->GetLowerRefYaxis()->SetRangeUser(0.5,1.5);
        rp->GetLowerRefYaxis()->SetTitle("Ratio(++/--)");

        TLegend* l =  rp->GetUpperPad()-> BuildLegend();
        l->Clear();
        l->AddEntry(h1,"mu+mu+","lep");
        l->AddEntry(h2,"mu-mu-","lep");

        rp->GetUpperPad()->Modified();
        rp->GetUpperPad()->Update();
        c1->Update();       
        c1->SaveAs(c);
        
        h1 = NULL;
        h2 = NULL;
        c1->Clear();
      }
    }
  }
  c1->Close();
}
