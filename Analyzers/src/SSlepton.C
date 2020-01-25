#include "SSlepton.h"

SSlepton::SSlepton(){

}

void SSlepton::initializeAnalyzer(){

  RunNI = HasFlag("RunNI");

  //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/SSlepton.h
  MuonIDs = {    
    "POGLoose",                            //DataFormat/src/Muon.C
    "POGMedium",
    "POGTight",
    "POGHighPt",
  };
  //==== corresponding Muon ID SF Keys for mcCorr->MuonID_SF()
  MuonIDSFKeys = {      
    "",                     //SKFlatAnalyzer/data/Run2Legacy_v3/2016/ID/histmap.txt
    "NUM_MediumID_DEN_genTracks",
    "NUM_TightID_DEN_genTracks",
    "NUM_HighPtID_DEN_genTracks",
  };

  //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
  //==== You can define sample-dependent or year-dependent variables here
  //==== (Example) Year-dependent variables
  //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/SSlepton.h 
  //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(Dataer==~~)" for every event (let's save cpu time).
  //==== Then, do it here, which only ran once for each macro

  if(DataYear==2016){
    IsoMuTriggerName = "HLT_IsoMu24_v";        //SKFlatAnalyzer/script/PDandTrigger
    TriggerSafePtCut = 26.;
  }
  else if(DataYear==2017){
    IsoMuTriggerName = "HLT_IsoMu27_v";
    TriggerSafePtCut = 29.;
  }

  else if(DataYear==2018){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }

  cout << "[SSlepton::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[SSlepton::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== Test btagging code
  //==== add taggers and WP that you want to use in analysis
  std::vector<Jet::Tagger> vtaggers;
  vtaggers.push_back(Jet::DeepCSV);

  std::vector<Jet::WP> v_wps;
  v_wps.push_back(Jet::Medium);

  //=== list of taggers, WP, setup systematics, use period SFs
  SetupBTagger(vtaggers,v_wps, true, true);

}

SSlepton::~SSlepton(){

  //==== Destructor of this Analyzer

}

void SSlepton::executeEvent(){

  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  AllJets = GetAllJets();

  weight_Prefire = GetPrefireWeight(0);

  AnalyzerParameter param;

  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

    TString MuonID = MuonIDs.at(it_MuonID);
    TString MuonIDSFKey = MuonIDSFKeys.at(it_MuonID);

    param.Clear();

    param.syst_ = AnalyzerParameter::Central;

    param.Name = MuonID;
    param.Muon_Tight_ID = MuonID;
    param.Muon_ID_SF_Key = MuonIDSFKey;
    param.Electron_Veto_ID = "passVetoID";
    param.Jet_ID = "tight";

    executeEventFromParameter(param);
  }
}

void SSlepton::executeEventFromParameter(AnalyzerParameter param){

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

  //==============
  //==== Trigger
  //==============
  if(! (ev.PassTrigger(IsoMuTriggerName) )) return;

  //======================
  //==== Copy AllObjects
  //======================

  vector<Muon> this_AllMuons = AllMuons;
  vector<Electron> this_AllElectrons = AllElectrons;
  vector<Jet> this_AllJets = AllJets;

  //=================================================
  //==== Then, apply ID selections using this_AllXXX
  //==================================================

  vector<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 20., 2.4);
  vector<Electron> eles = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 35., 2.5);

  vector<Jet> alljets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
  
  //=======================
  //==== Sort in pt-order
  //=======================

  //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(eles.begin(), eles.end(), PtComparing);
  //==== 2) alljets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
  std::sort(alljets.begin(), alljets.end(), PtComparing);

 
  //=========================
  //==== Event selections..
  //=========================

  //==== dimuon
  if(muons.size() != 2) return;

  //==== leading muon has trigger-safe pt
  if(muons.at(0).Pt() <= TriggerSafePtCut ) return;
 
  //==== 3rd lepton veto(only dimuon left)
  if(eles.size() != 0) return;

  //==== same sign dimuon 
  if(muons.at(0).Charge()*muons.at(1).Charge()<0) return;

  if(RunNI){
    if(!(muons.at(0).RelIso() < 0.15 && muons.at(1).RelIso() > 0.3) || (muons.at(0).RelIso() > 0.3 && muons.at(1).RelIso() < 0.15)) return;
  }
  else{
    if(!(muons.at(0).RelIso() < 0.15 && muons.at(1).RelIso() < 0.15)) return;
  }
  //===================
  //==== Event weight
  //===================

  double weight = 1.;
  //==== If MC
  if(!IsDATA){

    //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
    //==== So, you have to multiply trigger luminosity
    //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"

    weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");

    //==== MCweight is +1 or -1. Should be multiplied if you are using e.g., aMC@NLO NLO samples
    // steps that need to be done for each event. 
    weight *= ev.MCweight();

    //==== L1Prefire reweight
    weight *= weight_Prefire;

    //==== Example of applying Muon scale factors
    for(unsigned int i=0; i<muons.size(); i++){

      double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), weight, muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), weight, muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      weight *= this_idsf*this_isosf;

    }
  }  
  
  Charge_Plus(ev, param, weight, muons, eles, alljets);
  Charge_Minus(ev, param, weight, muons, eles, alljets);
  //Ratio(ev, param, weight, muons, eles, alljets);

}

void SSlepton::Charge_Plus(Event ev, AnalyzerParameter param, double weight,std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> alljets){

  TString sign = "plus";
  TString dir = param.Name + "/" + sign;

  Particle METv = ev.GetMETVector();
  double MET = METv.Pt();

  Particle ll  = muons.at(0) + muons.at(1);

  double mu0_pt, mu1_pt;
  mu0_pt = muons.at(0).Pt();
  mu1_pt = muons.at(1).Pt();

  vector<Jet> jets = JetsVetoLeptonInside(alljets, eles, muons);

  int Nbjet=0;

  for(unsigned int ij = 0 ; ij < alljets.size(); ij++){
    if(IsBTagged(alljets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) Nbjet++; // method for getting btag with SF applied to MC
  }
  
  if (muons.at(0).Charge() <  0)  return;

  FillHist(dir+"/mll", ll.M(), weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_pt", mu0_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu1_pt", mu1_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_eta", muons.at(0).Eta(), weight, 60, -3., 3.);
  FillHist(dir+"/mu1_eta", muons.at(1).Eta(), weight, 60, -3., 3.);

  if (MET < 40.) return;

  FillHist(dir+"/Njet", jets.size(), weight, 10, 0., 10.);
  FillHist(dir+"/Nbjet", Nbjet, weight, 10, 0.,10.);
   
  FillHist(dir+"/mll_MET40", ll.M(), weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_pt_MET40", mu0_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu1_pt_MET40", mu1_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_eta_MET40", muons.at(0).Eta(), weight, 60, -3., 3.);
  FillHist(dir+"/mu1_eta_MET40", muons.at(1).Eta(), weight, 60, -3., 3.);
 
  if (Nbjet == 0){
    FillHist(dir+"/mll_b_veto", ll.M(), weight, 300, 0., 3000.);
    FillHist(dir+"/mu0_pt_b_veto", mu0_pt, weight, 300, 0., 3000.);
    FillHist(dir+"/mu1_pt_b_veto", mu1_pt, weight, 300, 0., 3000.);
    FillHist(dir+"/mu0_eta_b_veto", muons.at(0).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mu1_eta_b_veto", muons.at(1).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mt1_b_veto", MT(muons.at(0), METv), weight, 50, 0., 500.);
    FillHist(dir+"/mt2_b_veto", MT(muons.at(1), METv), weight, 50, 0., 500.);  
  
    if (ll.M() < 40.){
      FillHist(dir+"/mu0_pt_b_veto_gg", mu0_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mu1_pt_b_veto_gg", mu1_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mt1_b_veto_gg", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_b_veto_gg", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
    else{
      FillHist(dir+"/mu0_pt_b_veto_W", mu0_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mu1_pt_b_veto_W", mu1_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mt1_b_veto_W", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_b_veto_W", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
  }

  else{
    FillHist(dir+"/Njet_1b", alljets.size(), weight, 10, 0., 10.);   
    FillHist(dir+"/mll_1b", ll.M(), weight, 3000, 0., 3000.);
    FillHist(dir+"/mu0_pt_1b", mu0_pt, weight, 3000, 0., 3000.);
    FillHist(dir+"/mu1_pt_1b", mu1_pt, weight, 3000, 0., 3000.);
    FillHist(dir+"/mu0_eta_1b", muons.at(0).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mu1_eta_1b", muons.at(1).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mt1_1b", MT(muons.at(0), METv), weight, 50, 0., 500.);
    FillHist(dir+"/mt2_1b", MT(muons.at(1), METv), weight, 50, 0., 500.);
 
    if (ll.M() < 40.){
      FillHist(dir+"/mu0_pt_1b_gg", mu0_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mu1_pt_1b_gg", mu1_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mt1_1b_gg", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_1b_gg", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
    else{
      FillHist(dir+"/mu0_pt_1b_W", mu0_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mu1_pt_1b_W", mu1_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mt1_1b_W", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_1b_W", MT(muons.at(1), METv), weight, 50, 0., 500.);

    }
  }
}

void SSlepton::Charge_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> alljets){

  TString sign = "minus";
  TString dir = param.Name + "/" + sign;

  Particle METv = ev.GetMETVector();
  double MET = METv.Pt();

  Particle ll  = muons.at(0) + muons.at(1);

  double mu0_pt, mu1_pt;
  mu0_pt = muons.at(0).Pt();
  mu1_pt = muons.at(1).Pt();

  vector<Jet> jets = JetsVetoLeptonInside(alljets, eles, muons);

  int Nbjet=0;
  for(unsigned int ij = 0 ; ij < alljets.size(); ij++){
    if(IsBTagged(alljets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) Nbjet++; // method for getting btag with SF applied to MC
  }

  if (muons.at(0).Charge() > 0) return;

  FillHist(dir+"/mll", ll.M(), weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_pt", mu0_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu1_pt", mu1_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_eta", muons.at(0).Eta(), weight, 60, -3., 3.);
  FillHist(dir+"/mu1_eta", muons.at(1).Eta(), weight, 60, -3., 3.);

  if (MET < 40.) return;

  FillHist(dir+"/Njet", jets.size(), weight, 10, 0., 10.);
  FillHist(dir+"/Nbjet", Nbjet, weight, 10, 0.,10.);
   
  FillHist(dir+"/mll_MET40", ll.M(), weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_pt_MET40", mu0_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu1_pt_MET40", mu1_pt, weight, 300, 0., 3000.);
  FillHist(dir+"/mu0_eta_MET40", muons.at(0).Eta(), weight, 60, -3., 3.);
  FillHist(dir+"/mu1_eta_MET40", muons.at(1).Eta(), weight, 60, -3., 3.);
 
  if (Nbjet == 0){
    FillHist(dir+"/mll_b_veto", ll.M(), weight, 300, 0., 3000.);
    FillHist(dir+"/mu0_pt_b_veto", mu0_pt, weight, 300, 0., 3000.);
    FillHist(dir+"/mu1_pt_b_veto", mu1_pt, weight, 300, 0., 3000.);
    FillHist(dir+"/mu0_eta_b_veto", muons.at(0).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mu1_eta_b_veto", muons.at(1).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mt1_b_veto", MT(muons.at(0), METv), weight, 50, 0., 500.);
    FillHist(dir+"/mt2_b_veto", MT(muons.at(1), METv), weight, 50, 0., 500.);  
  
    if (ll.M() < 40.){
      FillHist(dir+"/mu0_pt_b_veto_gg", mu0_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mu1_pt_b_veto_gg", mu1_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mt1_b_veto_gg", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_b_veto_gg", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
    else{
      FillHist(dir+"/mu0_pt_b_veto_W", mu0_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mu1_pt_b_veto_W", mu1_pt, weight, 300, 0., 3000.);
      FillHist(dir+"/mt1_b_veto_W", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_b_veto_W", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
  }

  else{
    FillHist(dir+"/Njet_1b", alljets.size(), weight, 10, 0., 10.);   
    FillHist(dir+"/mll_1b", ll.M(), weight, 3000, 0., 3000.);
    FillHist(dir+"/mu0_pt_1b", mu0_pt, weight, 3000, 0., 3000.);
    FillHist(dir+"/mu1_pt_1b", mu1_pt, weight, 3000, 0., 3000.);
    FillHist(dir+"/mu0_eta_1b", muons.at(0).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mu1_eta_1b", muons.at(1).Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/mt1_1b", MT(muons.at(0), METv), weight, 50, 0., 500.);
    FillHist(dir+"/mt2_1b", MT(muons.at(1), METv), weight, 50, 0., 500.);
 
    if (ll.M() < 40.){
      FillHist(dir+"/mu0_pt_1b_gg", mu0_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mu1_pt_1b_gg", mu1_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mt1_1b_gg", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_1b_gg", MT(muons.at(1), METv), weight, 50, 0., 500.);
    }
    else{
      FillHist(dir+"/mu0_pt_1b_W", mu0_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mu1_pt_1b_W", mu1_pt, weight, 3000, 0., 3000.);
      FillHist(dir+"/mt1_1b_W", MT(muons.at(0), METv), weight, 50, 0., 500.);
      FillHist(dir+"/mt2_1b_W", MT(muons.at(1), METv), weight, 50, 0., 500.);

    }
  }
}
