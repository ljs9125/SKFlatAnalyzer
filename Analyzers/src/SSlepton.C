#include "SSlepton.h"

SSlepton::SSlepton(){

}

void SSlepton::initializeAnalyzer(){

  RunNI = HasFlag("RunNI");

  //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/SSlepton.h
  MuonIDs = {    
    "POGLoose",                            //DataFormat/src/Muon.C
    "POGTight",
  };
  //==== corresponding Muon ID SF Keys for mcCorr->MuonID_SF()
  MuonIDSFKeys = {      
    "",                     //SKFlatAnalyzer/data/Run2Legacy_v3/2016/ID/histmap.txt
    "NUM_TightID_DEN_genTracks",
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

  //================================================================
  //====  Example 1
  //====  llon Z-peak events with two muon IDs, with systematics
  //================================================================

  //==== *IMPORTANT TO SAVE CPU TIME*
  //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
  //==== and then check ID booleans.
  //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
  //==== We are now running systematics, and you don't want to do this for every systematic sources
  //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/SSlepton.h,
  //==== and save muons objects at the very beginning of executeEvent().
  //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
  AllMuons = GetAllMuons();
  AllElectrons = GetAllElectrons();
  //=== Jets too
  AllJets = GetAllJets();

  //==== Get L1Prefire reweight
  //==== If data, 1.;
  //==== If MC && DataYear > 2017, 1.;
  //==== If MC && DataYear <= 2017, we have to reweight the event with this value
  //==== I defined "double weight_Prefire;" in Analyzers/include/SSlepton.h
  weight_Prefire = GetPrefireWeight(0);

  //==== Declare AnalyzerParameter

  AnalyzerParameter param;

  //==== Loop over muon IDs

  for(unsigned int it_MuonID=0; it_MuonID<MuonIDs.size(); it_MuonID++){

    TString MuonID = MuonIDs.at(it_MuonID);
    TString MuonIDSFKey = MuonIDSFKeys.at(it_MuonID);

    //==== 1) First, let's run Central values of the systematics

    //==== clear parameter set
    param.Clear();

    //==== set which systematic sources you want to run this time
    //==== default syst_ is AnalyzerParameter::Central
    param.syst_ = AnalyzerParameter::Central;

    //==== set name of the parameter set
    //==== this will be used for the directory name of histograms
    param.Name = MuonID;

    //==== You can define lepton ID string here
    param.Muon_Tight_ID = MuonID;
    param.Muon_ID_SF_Key = MuonIDSFKey;
    param.Electron_Veto_ID = "passVetoID";

    //==== And, Jet ID
    param.Jet_ID = "tight";

    //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set
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
  vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);

  //=======================
  //==== Sort in pt-order
  //=======================

  //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(eles.begin(), eles.end(), PtComparing);
  //==== 2) jets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
  std::sort(jets.begin(), jets.end(), PtComparing);

 
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

      double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      weight *= this_idsf*this_isosf;

    }
  }  
  
  Charge_Plus(ev, param, weight, muons, eles, jets);
  Charge_Minus(ev, param, weight, muons, eles, jets);

}

void SSlepton::Charge_Plus(Event ev, AnalyzerParameter param, double weight,std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets){

  Particle METv = ev.GetMETVector();
  double MET = METv.Pt();

  Particle ll  = muons.at(0) + muons.at(1);

  double Mu_1st_pt, Mu_2nd_pt;
  Mu_1st_pt = muons.at(0).Pt();
  Mu_2nd_pt = muons.at(1).Pt();

  int Nbjet=0;

  for(unsigned int ij = 0 ; ij < jets.size(); ij++){
    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) Nbjet++; // method for getting btag with SF applied to MC
  }
  
  if (muons.at(0).Charge() <  0)  return;

  FillHist(param.Name+"/plus/dimu_invmass", ll.M(), weight, 3000, 0., 3000.);
  FillHist(param.Name+"/plus/mu_1st_pt", Mu_1st_pt, weight, 3000, 0., 3000.);
  FillHist(param.Name+"/plus/mu_2nd_pt", Mu_2nd_pt, weight, 3000, 0., 3000.);

  if (MET < 40.) return;

  FillHist(param.Name+"/plus/Njet", jets.size(), weight, 10, 0., 10.);
  FillHist(param.Name+"/plus/Nbjet", Nbjet, weight, 10, 0.,10.);
   
  FillHist(param.Name+"/plus/dimu_invmass_MET40", ll.M(), weight, 3000, 0., 3000.);
  FillHist(param.Name+"/plus/mu_1st_pt_MET40", Mu_1st_pt, weight, 3000, 0., 3000.);
  FillHist(param.Name+"/plus/mu_2nd_pt_MET40", Mu_2nd_pt, weight, 3000, 0., 3000.);

  if (Nbjet == 0){
    FillHist(param.Name+"/plus/dimu_invmass_Nobjet", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_1st_pt_Nobjet", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_2nd_pt_Nobjet", Mu_2nd_pt, weight, 3000, 0., 3000.);
  }
  else{
    FillHist(param.Name+"/plus/Njet_1b", jets.size(), weight, 10, 0., 10.);   

    FillHist(param.Name+"/plus/dimu_invmass_1b", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_1st_pt_1b", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_2nd_pt_1b", Mu_2nd_pt, weight, 3000, 0., 3000.);

    if(jets.size() == 0) return;

    FillHist(param.Name+"/plus/dimu_invmass_1b1j", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_1st_pt_1b1j", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/plus/mu_2nd_pt_1b1j", Mu_2nd_pt, weight, 3000, 0., 3000.); 
  }
}

void SSlepton::Charge_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets){

  Particle METv = ev.GetMETVector();
  double MET = METv.Pt();

  Particle ll  = muons.at(0) + muons.at(1);

  double Mu_1st_pt, Mu_2nd_pt;
  Mu_1st_pt = muons.at(0).Pt();
  Mu_2nd_pt = muons.at(1).Pt();

  int Nbjet=0;
  for(unsigned int ij = 0 ; ij < jets.size(); ij++){
    if(IsBTagged(jets.at(ij), Jet::DeepCSV, Jet::Medium,true,0)) Nbjet++; // method for getting btag with SF applied to MC
  }

  if (muons.at(0).Charge() > 0) return;

  FillHist(param.Name+"/minus/dimu_invmass", ll.M(), weight, 3000, 0., 3000.);
  FillHist(param.Name+"/minus/mu_1st_pt", Mu_1st_pt, weight, 3000, 0., 3000.);
  FillHist(param.Name+"/minus/mu_2nd_pt", Mu_2nd_pt, weight, 3000, 0., 3000.);

  if (MET < 40.) return;

  FillHist(param.Name+"/minus/Njet", jets.size(), weight, 10, 0., 10.);
  FillHist(param.Name+"/minus/Nbjet", Nbjet, weight, 10, 0., 10.);  

  FillHist(param.Name+"/minus/dimu_invmass_MET40", ll.M(), weight, 3000, 0., 3000.);
  FillHist(param.Name+"/minus/mu_1st_pt_MET40", Mu_1st_pt, weight, 3000, 0., 3000.);
  FillHist(param.Name+"/minus/mu_2nd_pt_MET40", Mu_2nd_pt, weight, 3000, 0., 3000.);

  if (Nbjet == 0){
    FillHist(param.Name+"/minus/dimu_invmass_Nobjet", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_1st_pt_Nobjet", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_2nd_pt_Nobjet", Mu_2nd_pt, weight, 3000, 0., 3000.);
  }
  else{
    FillHist(param.Name+"/minus/Njet", jets.size(), weight, 10, 0., 10.);

    FillHist(param.Name+"/minus/dimu_invmass_1b", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_1st_pt_1b", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_2nd_pt_1b", Mu_2nd_pt, weight, 3000, 0., 3000.);

    if(jets.size() == 0) return;
    FillHist(param.Name+"/minus/dimu_invmass_1b1j", ll.M(), weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_1st_pt_1b1j", Mu_1st_pt, weight, 3000, 0., 3000.);
    FillHist(param.Name+"/minus/mu_2nd_pt_1b1j", Mu_2nd_pt, weight, 3000, 0., 3000.);
  }
}



