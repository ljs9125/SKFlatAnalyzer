#ifndef SSlepton_h
#define SSlepton_h

#include "AnalyzerCore.h"

class SSlepton : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunNI;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Electron> AllElectrons;
  vector<Jet> AllJets;

  double weight_Prefire;

  SSlepton();
  ~SSlepton();

  void Iso_Plus(Event ev, AnalyzerParameter param, double weight, std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets);
  void Iso_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets); 
  void NIso_Plus(Event ev, AnalyzerParameter param, double weight, std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets);
  void NIso_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets); 
  void NNIso_Plus(Event ev, AnalyzerParameter param, double weight, std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets);
  void NNIso_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets); 


  void Plot_All(TString dir, std::vector<Muon> muons, Particle ll, Particle METv, std::vector<Jet> jets, std::vector<Jet> alljets, std::vector<Jet> bjet, int Nbjet, double weight);
  void FillMuonPlots(vector<Muon> muons, TString this_dir, TString this_region , double weight);
  void FillJetsPlots(vector<Jet> jets, vector<Jet> bjet, TString this_dir, TString this_region, double weight);
};



#endif

