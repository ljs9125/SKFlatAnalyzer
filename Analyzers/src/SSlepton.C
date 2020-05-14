#include "SSlepton.h"

SSlepton::SSlepton(){

}

void SSlepton::initializeAnalyzer(){

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
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }

  else if(DataYear==2018){
    IsoMuTriggerName = "HLT_IsoMu24_v";
    TriggerSafePtCut = 26.;
  }

  cout << "[SSlepton::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
  cout << "[SSlepton::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;

  //==== B-Tagging
  //==== add taggers and WP that you want to use in analysis
  std::vector<JetTagging::Parameters> jtps;
  //==== If you want to use 1a or 2a method,
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
  //==== set
  mcCorr->SetJetTaggingParameters(jtps);

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
    param.Muon_Veto_ID = "POGLoose";
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
  vector<Muon> muons_veto = SelectMuons(this_AllMuons, param.Muon_Veto_ID, 20., 2.4);
  vector<Electron> eles = SelectElectrons(this_AllElectrons, param.Electron_Veto_ID, 35., 2.5);

  vector<Jet> alljets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
  
  //=======================
  //==== Sort in pt-order
  //=======================

  //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
  std::sort(muons.begin(), muons.end(), PtComparing);
  std::sort(muons_veto.begin(), muons_veto.end(), PtComparing);
  std::sort(eles.begin(), eles.end(), PtComparing);
  //==== 2) alljets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
  std::sort(alljets.begin(), alljets.end(), PtComparing);

 
  //=========================
  //==== Event selections..
  //=========================

  //==== dimuon
  if(muons.size() != 2) return;

  //==== 3rd lepton veto(only dimuon)
  if(eles.size() != 0) return;
  if(muons.size() != muons_veto.size()) return;

  //==== leading muon has trigger-safe pt
  if(muons.at(0).Pt() <= TriggerSafePtCut ) return; 

  //==== same sign dimuon 
  if(muons.at(0).Charge()*muons.at(1).Charge()<0) return;

  vector<Jet> bjet_L;
  vector<Jet> bjet_M;
  vector<Jet> bjet_T;

  int Nbjet_L=0;
  int Nbjet_M=0;
  int Nbjet_T=0;

  for(unsigned int ij = 0 ; ij < alljets.size(); ij++){

    double this_discr = alljets.at(ij).GetTaggerResult(JetTagging::DeepCSV);

    if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)){
      Nbjet_M++;
      bjet_M.push_back(alljets.at(ij));
    }
    if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Loose)){
      Nbjet_L++;
      bjet_L.push_back(alljets.at(ij));
    }
    if( this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Tight)){
      Nbjet_T++;
      bjet_T.push_back(alljets.at(ij));
    }
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
    /*
    for(unsigned int i=0; i<muons.size(); i++){

      double this_idsf  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(i).Eta(), weight, muons.at(i).MiniAODPt());

      //==== If you have iso SF, do below. Here we don't.
      //double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(i).Eta(), weight, muons.at(i).MiniAODPt());
      double this_isosf = 1.;

      weight *= this_idsf*this_isosf;

    } */
  }  

  TString iso;
  TString sign;
  TString dir = param.Name + "/";
  TString temp_dir; 

  Particle METv = ev.GetMETVector();
  double MET = METv.Pt();

  Particle ll  = muons.at(0) + muons.at(1);

  //==== reject collinear div
  if (ll.M() < 10) return;
  //==== reduce QCD bkgd.
  if (MET < 40.) return;

   //==== assign charge 
  if (muons.at(0).Charge() < 0) sign = "minus";
  else if (muons.at(0).Charge() > 0) sign = "plus";
  else return;

  //==== Divide regions according to Muon isolation
  if(muons.at(0).RelIso() < 0.15 && muons.at(1).RelIso() < 0.15) iso = "II";
  else if((muons.at(0).RelIso() < 0.15 && muons.at(1).RelIso() > 0.3) || (muons.at(0).RelIso() > 0.3 && muons.at(1).RelIso() < 0.15)) iso = "NI";
  else if((muons.at(0).RelIso() > 0.3 && muons.at(1).RelIso() > 0.3)) iso = "NN";
  else return;

  dir = dir + iso + "/" + sign;
  temp_dir = dir;

  cout << "current directory: " << dir << endl;

  dir = dir + "/SortMuByPt/Loose";
  Plot_All(dir, muons, ll, METv, alljets, bjet_L, Nbjet_L, weight);
  dir = temp_dir + "/SortMuByPt/Medium";
  Plot_All(dir, muons, ll, METv, alljets, bjet_M, Nbjet_M, weight);
  dir = temp_dir + "/SortMuByPt/Tight";
  Plot_All(dir, muons, ll, METv, alljets, bjet_T, Nbjet_T, weight);


  dir = temp_dir;

  if(dir.Contains("NI")){
    dir = dir + "/SortMuByIso";
    temp_dir = dir;

    vector<Muon> IsoMu;
    if(muons.at(0).RelIso() < 0.15){
      IsoMu.push_back(muons.at(0));
      IsoMu.push_back(muons.at(1));
    }  
    else if(muons.at(0).RelIso() > 0.3){
      IsoMu.push_back(muons.at(1));
      IsoMu.push_back(muons.at(0));
    }

    dir = dir + "/Loose";
    Plot_All(dir, IsoMu, ll, METv, alljets, bjet_L, Nbjet_L, weight);
    dir = temp_dir + "/Medium";
    Plot_All(dir, IsoMu, ll, METv, alljets, bjet_M, Nbjet_M, weight);
    dir = temp_dir + "/Tight";
    Plot_All(dir, IsoMu, ll, METv, alljets, bjet_T, Nbjet_T, weight);

    
  }

}


void SSlepton::FillMuonPlots(std::vector<Muon> muons, TString this_dir, TString this_region, double weight){
  
  for(unsigned int i=0; i<muons.size(); i++){

    TString this_itoa = TString::Itoa(i,10);

    Muon muon = muons[i];

    FillHist(this_dir+"/Muon_"+this_itoa+"_Pt_"+this_region, muon.Pt(), weight, 1000, 0., 1000.);
    FillHist(this_dir+"/Muon_"+this_itoa+"_Eta_"+this_region, muon.Eta(), weight, 60, -3., 3.);
    FillHist(this_dir+"/Muon_"+this_itoa+"_RelIso_"+this_region, muon.RelIso(), weight, 100, 0., 1.);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_MiniRelIso_"+this_region, muon.MiniRelIso(), weight, 100, 0., 1.);

    //FillHist(this_dir+"/Muon_"+this_itoa+"_dXY_"+this_region, fabs(muon.dXY()), weight, 500, 0., 0.05);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_dXYSig_"+this_region, fabs(muon.dXY()/muon.dXYerr()), weight, 100, 0., 10);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_dZ_"+this_region, fabs(muon.dZ()), weight, 500, 0., 0.5);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_dZSig_"+this_region, fabs(muon.dZ()/muon.dZerr()), weight, 100, 0., 10);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_IP3D_"+this_region, fabs(muon.IP3D()), weight, 500, 0., 0.5);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_IP3DSig_"+this_region, fabs(muon.IP3D()/muon.IP3Derr()), weight, 100, 0., 10);
    //FillHist(this_dir+"/Muon_"+this_itoa+"_Chi2_"+this_region, muon.Chi2(), weight, 500, 0., 50.);
    FillHist(this_dir+"/Muon_"+this_itoa+"_TrkRelIso_"+this_region, muon.TrkIso()/muon.TuneP4().Pt(), weight, 100, 0., 1.);
    
  }
}

void SSlepton::FillJetsPlots(std::vector<Jet> jets, std::vector<Jet> bjet, TString this_dir, TString this_region, double weight){

  for(unsigned int i=0; i<jets.size(); i++){

    TString this_itoa = TString::Itoa(i,10);
    FillHist(this_dir+"/Jet_"+this_itoa+"_Pt_"+this_region, jets.at(i).Pt(), weight, 1000, 0., 1000.);
    FillHist(this_dir+"/Jet_"+this_itoa+"_Eta_"+this_region, jets.at(i).Eta(), weight, 60, -3., 3.);

  }

  for(unsigned int i=0; i<bjet.size(); i++){

    TString this_itoa = TString::Itoa(i,10);
    FillHist(this_dir+ "/Bjet_"+this_itoa+"_Pt_"+this_region, bjet.at(i).Pt(), weight, 1000, 0., 1000.);
    FillHist(this_dir+"/Bjet_"+this_itoa+"_Eta_"+this_region, bjet.at(i).Eta(), weight, 60, -3., 3.);

  }
}

void SSlepton::Plot_All(TString dir, std::vector<Muon> muons, Particle ll, Particle METv, std::vector<Jet> alljets, std::vector<Jet> bjet , int Nbjet ,double weight){

  if (Nbjet == 0){
    dir = dir + "/b_veto";
    FillHist(dir+"/Njet_b_veto", alljets.size(), weight, 10, 0., 10.);
    FillHist(dir+"/mll_b_veto", ll.M(), weight, 3000, 0., 3000.);
    FillMuonPlots(muons, dir,"b_veto", weight);
    FillHist(dir+"/METv_pt_b_veto", METv.Pt(), weight, 1000, 0., 1000.);
    FillHist(dir+"/METv_eta_b_veto", METv.Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/METv_phi_b_veto", METv.Phi(), weight, 100, -5., 5.);
    FillHist(dir+"/mt1_b_veto", MT(muons.at(0), METv), weight, 500, 0., 500.);
    FillHist(dir+"/mt2_b_veto", MT(muons.at(1), METv), weight, 500, 0., 500.);  
    FillHist(dir+"/dphi1_METv_b_veto", muons.at(0).DeltaPhi(METv), weight, 100, -5., 5.);
    FillHist(dir+"/dphi2_METv_b_veto", muons.at(1).DeltaPhi(METv), weight, 100, -5., 5.);
    FillJetsPlots(alljets, bjet, dir, "b_veto", weight);

  }
  else{
    dir = dir + "/1b";
    FillHist(dir+"/Njet_1b", alljets.size(), weight, 10, 0., 10.);
    FillHist(dir+"/mll_1b", ll.M(), weight, 3000, 0., 3000.);
    FillMuonPlots(muons, dir, "1b", weight);
    FillHist(dir+"/METv_pt_1b", METv.Pt(), weight, 1000, 0., 1000.);
    FillHist(dir+"/METv_eta_1b", METv.Eta(), weight, 60, -3., 3.);
    FillHist(dir+"/METv_phi_1b", METv.Phi(), weight, 100, -5., 5.);
    FillHist(dir+"/mt1_1b", MT(muons.at(0), METv), weight, 500, 0., 500.);
    FillHist(dir+"/mt2_1b", MT(muons.at(1), METv), weight, 500, 0., 500.);
    FillHist(dir+"/dphi1_METv_1b", muons.at(0).DeltaPhi(METv), weight, 100, -5., 5.);
    FillHist(dir+"/dphi2_METv_1b", muons.at(1).DeltaPhi(METv), weight, 100, -5., 5.);
    FillJetsPlots(alljets, bjet, dir, "1b", weight);
    FillHist(dir+"/deltaR_b_mu1_1b", muons.at(0).DeltaR(bjet.at(0)), weight, 1000, 0., 10.);
    FillHist(dir+"/deltaR_b_mu2_1b", muons.at(1).DeltaR(bjet.at(0)), weight, 1000, 0., 10.);

    if (muons.at(1).DeltaR(bjet.at(0)) < 0.3) {
      dir = dir + "/MuonInsideBjet" ;
      FillMuonPlots(muons, dir, "1b", weight);
      FillHist(dir+"/METv_pt_1b", METv.Pt(), weight, 1000, 0., 1000.);
      FillHist(dir+"/METv_eta_1b", METv.Eta(), weight, 60, -3., 3.);
      FillHist(dir+"/METv_phi_1b", METv.Phi(), weight, 100, -5., 5.);
      FillHist(dir+"/mt1_1b", MT(muons.at(0), METv), weight, 500, 0., 500.);
      FillHist(dir+"/mt2_1b", MT(muons.at(1), METv), weight, 500, 0., 500.);
      FillHist(dir+"/dphi1_METv_1b", muons.at(0).DeltaPhi(METv), weight, 100, -5., 5.);
      FillHist(dir+"/dphi2_METv_1b", muons.at(1).DeltaPhi(METv), weight, 100, -5., 5.);
      if ( MT(muons.at(1), METv) > 50 && MT(muons.at(1), METv) < 110) {
        dir = dir + "/MT_GT_50_LT_110" ;
        Particle W_rel = muons.at(1) + METv;
        FillHist(dir+"/dphi_NonIso_METv_1b", W_rel.DeltaPhi(bjet.at(0)), weight, 100, -5., 5.);
        FillHist(dir+"/deltaR_NonIso_METv_1b", W_rel.DeltaR(bjet.at(0)), weight, 1000, 0., 10.);
      }
    }

  }

/*  dir = temp;

  if (Nbjet_L == 0){
    dir = dir + "/Loose_b_veto";
    FillHist(dir+"/mll_Loose_b_veto", ll.M(), weight, 3000, 0., 3000.);
    FillHist(dir+"/mt1_Loose_b_veto", MT(muons.at(0), METv), weight, 500, 0., 500.);
    FillHist(dir+"/mt2_Loose_b_veto", MT(muons.at(1), METv), weight, 500, 0., 500.);
    FillHist(dir+"/dphi1_Loose_b_veto", muons.at(0).DeltaPhi(METv), weight, 100, -5., 5.);
    FillHist(dir+"/dphi2_Loose_b_veto", muons.at(1).DeltaPhi(METv), weight, 100, -5., 5.);

  }
  if (Nbjet_T > 0){
    dir = dir + "/Tight_1b";
    FillHist(dir+"/mll_Tight_1b", ll.M(), weight, 3000, 0., 3000.);
    FillHist(dir+"/mt1_Tight_1b", MT(muons.at(0), METv), weight, 500, 0., 500.);
    FillHist(dir+"/mt2_Tight_1b", MT(muons.at(1), METv), weight, 500, 0., 500.);
    FillHist(dir+"/dphi1_METv_Tight_1b", muons.at(0).DeltaPhi(METv), weight, 100, -5., 5.);
    FillHist(dir+"/dphi2_METv_Tight_1b", muons.at(1).DeltaPhi(METv), weight, 100, -5., 5.);
    FillHist(dir+"/deltaR_b_mu1_1b", muons.at(0).DeltaR(bjet_T.at(0)), weight, 1000, 0., 10.);
    FillHist(dir+"/deltaR_b_mu2_1b", muons.at(1).DeltaR(bjet_T.at(0)), weight, 1000, 0., 10.);
  }
*/

}
