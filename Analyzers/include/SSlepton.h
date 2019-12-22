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

  void Charge_Plus(Event ev, AnalyzerParameter param, double weight, std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets);
  void Charge_Minus(Event ev, AnalyzerParameter param, double weight,  std::vector<Muon> muons, std::vector<Electron> eles, std::vector<Jet> jets);
};



#endif

