void SS_dimuon() {

  gStyle->SetOptLogy(1); 
  //gStyle->SetOptLogx(1); 
 
  TFile *DATA = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/DATA/SSlepton_SingleMuon.root");

  TFile *WZ_pythia = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_WZ_pythia.root");
  TFile *ZZ_pythia = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_ZZ_pythia.root");
 // TFile *WW_pythia = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_WW_pythia.root");


  TFile *DYJets = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_DYJets.root");

  TFile *TTLL_powheg = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_TTLL_powheg.root");
  TFile *TTLJ_powheg = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_TTLJ_powheg.root");
  TFile *WJets_MG = TFile::Open("/home/snuintern1/skflat/SKFlatOutput/Run2Legacy_v3/SSlepton/2016/SSlepton_WJets_MG.root");


  TCanvas *c = new TCanvas("c","c");

  TH1F *h1 = (TH1F*)DATA->Get("POGMedium_Central/SS_DATA_POGMedium_Central");

  TH1F *h2 = (TH1F*)WZ_pythia->Get("POGMedium_Central/SS_WZ_pythia_POGMedium_Central");
  TH1F *h3 = (TH1F*)ZZ_pythia->Get("POGMedium_Central/SS_ZZ_pythia_POGMedium_Central");
 // TH1F *h4 = (TH1F*)WW_pythia->Get("POGMedium_Central/SS_WW_pythia_POGMedium_Central");

  TH1F *h5 = (TH1F*)DYJets->Get("POGMedium_Central/SS_DYJets_POGMedium_Central");
 
  TH1F *h6 = (TH1F*)TTLJ_powheg->Get("POGMedium_Central/SS_TTLJ_powheg_POGMedium_Central");
  TH1F *h7 = (TH1F*)TTLL_powheg->Get("POGMedium_Central/SS_TTLL_powheg_POGMedium_Central");
  TH1F *h8 = (TH1F*)WJets_MG->Get("POGMedium_Central/SS_WJets_MG_POGMedium_Central");
  
  THStack *hs = new THStack("hs",""); 


  //draw data 
  c->cd();

  h1->Rebin(10);
  h1->Draw("E P");
  h1->GetXaxis()->SetRangeUser(10.,500.);
  h1->SetMarkerStyle(8);

  //draw MC Samples 1)charge mismatch 2)misidentification 3)prompt SS
  //Using THStack for MC Sample to compare with data
  
  h3->Add(h2);

  h3->Rebin(10);
  h3->Draw("HIST SAME");
  h3->SetLineColor(kGreen);
  h3->SetFillColor(kGreen);
  h3->GetXaxis()->SetRangeUser(10.,500.);
 
  //Use DY sample after Gen matching  

  h5->Rebin(10);
  h5->Draw("HIST SAME");
  h5->SetLineColor(kBlue);
  h5->SetFillColor(kBlue);
  
  h8->Add(h6);
  h8->Add(h7);

  h8->Rebin(10);
  h8->Draw("HIST SAME");
  h8->SetLineColor(kRed);
  h8->SetFillColor(kRed);

  hs->Add(h3);
  hs->Add(h5);
  hs->Add(h8);
  hs->Draw("SAME HIST");
 
//Axis Label
  h1->GetYaxis()->SetTitle("Entries");
  h1->GetXaxis()->SetTitle("M(mu+mu-) [GeV]");
  h1->SetTitle("SameSign dimuon invariant mass distribution"); 

  //Legend
  TLegend *l = new TLegend(0.2,0.3,0.4,0.4);
  l->AddEntry(h1,"DATA","lep"); 
  l->AddEntry(h3, "Prompt", "f");
  l->AddEntry(h5, "Non_prompt", "f");
  l->AddEntry(h8, "Misidentified", "f");
  l->Draw();
  
  //Save as pdf file
  c->SaveAs("./output/SSdimuon.pdf");

