void NonIso_2016() {


  //Set log scale 
  gStyle->SetOptLogy(1); 
 // gStyle->SetOptLogx(1); 
  
  //Open root file 
  TFile *DATA = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/NonIsoSec/2016/DATA/NonIsoSec_SingleMuon.root");

  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  TCanvas *c3 = new TCanvas("c3","c3");

  TCanvas *c4 = new TCanvas("c4","c4");
  TCanvas *c5 = new TCanvas("c5","c5");
  TCanvas *c6 = new TCanvas("c6","c6");

  TH1F *h1 = (TH1F*)DATA->Get("NonIso/POGLoose_Central/Mll_pp_DATA");
  TH1F *h2 = (TH1F*)DATA->Get("NonIsoPOGLoose_Central/Mll_mm_DATA");

  TH1F *h4 = (TH1F*)DATA->Get("NonIso/POGTight_Central/Mll_pp_DATA");
  TH1F *h5 = (TH1F*)DATA->Get("NonIsoPOGTight_Central/Mll_mm_DATA");

  //draw data 
  
  c1->cd();
  h1->Rebin(10);
  h1->Draw("E P");
  h1->GetXaxis()->SetRangeUser(10.,3000.);
  h1->SetMarkerStyle(8); 
 
  c2->cd();
  h2->Rebin(10);
  h2->Draw("E P");
  h2->GetXaxis()->SetRangeUser(10.,3000.);
  h2->SetMarkerStyle(8);

  c4->cd();
  h4->Rebin(10);
  h4->Draw("E P");
  h4->GetXaxis()->SetRangeUser(10.,3000.);
  h4->SetMarkerStyle(8);

  c5->cd();
  h5->Rebin(10);
  h5->Draw("E P");
  h5->GetXaxis()->SetRangeUser(10.,3000.);
  h5->SetMarkerStyle(8);

  //Becareful if you REBIN, you have to clone after Rebin
  TH1F *h3 = (TH1F*)h1-> Clone ("Ratio PP/MM"); 
  TH1F *h6 = (TH1F*)h4-> Clone ("Ratio PP/MM"); 


  c3->cd();
  c3->SetLogy(0);  //SetLog is function of TCanvas, SetOptLog for TStyle
  h3->Sumw2();
  h3->Divide(h2);
  h3->GetXaxis()->SetRangeUser(10.,1000.);
  h3->Draw("E P");
  h3->SetMarkerStyle(8);

  c6->cd();
  c6->SetLogy(0);  //SetLog is function of TCanvas, SetOptLog for TStyle
  h6->Sumw2();
  h6->Divide(h5);
  h6->GetXaxis()->SetRangeUser(10.,1000.);
  h6->Draw("E P");
  h6->SetMarkerStyle(8);



  //Axis Label

  h1->GetYaxis()->SetTitle("Entries");
  h1->GetXaxis()->SetTitle("M(mu+mu+) [GeV]");
  h1->SetTitle("++ dimuon invariant mass"); 

  h2->GetXaxis()->SetTitle("M(mu-mu-) [GeV]");
  h2->GetYaxis()->SetTitle("Entries");
  h2->SetTitle("-- invariant mass");
  
  h3->GetYaxis()->SetTitle("Entries");
  h3->GetXaxis()->SetTitle("M(ll) [GeV]");
  h3->SetTitle("Ratio plot ++/--"); 

  h4->GetYaxis()->SetTitle("Entries");
  h4->GetXaxis()->SetTitle("M(mu+mu+) [GeV]");
  h4->SetTitle("++ dimuon invariant mass"); 

  h5->GetXaxis()->SetTitle("M(mu-mu-) [GeV]");
  h5->GetYaxis()->SetTitle("Entries");
  h5->SetTitle("-- invariant mass");
  
  h6->GetYaxis()->SetTitle("Entries");
  h6->GetXaxis()->SetTitle("M(ll) [GeV]");
  h6->SetTitle("Ratio plot ++/--"); 



  //Fit ratio plot 
  h3->Fit("pol0","","",10.,1000.);
  c3->Update();

  h6->Fit("pol0","","",10.,1000.);

//Legend
//  TLegend *l = new TLegend(0.2,0.3,0.4,0.4);
//  l->AddEntry(h1,"DATA","lep"); 
//  l->Draw();
  
  //Save as pdf file
  c1->SaveAs("./output/2016/Loose/pp_dimuon_DATA.pdf");
  c2->SaveAs("./output/2016/Loose/mm_dimuon_DATA.pdf");
  c3->SaveAs("./output/2016/Loose/ratio_dimuon_DATA.pdf");
  c4->SaveAs("./output/2016/Tight/pp_dimuon_DATA.pdf");
  c5->SaveAs("./output/2016/Tight/mm_dimuon_DATA.pdf");
  c6->SaveAs("./output/2016/Tight/ratio_dimuon_DATA.pdf");
 
}
