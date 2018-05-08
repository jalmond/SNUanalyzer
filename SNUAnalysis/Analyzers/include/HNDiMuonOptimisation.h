#ifndef HNDiMuonOptimisation_h
#define HNDiMuonOptimisation_h

#include "AnalyzerCore.h"

class HNDiMuonOptimisation : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNDiMuonOptimisation();
  ~HNDiMuonOptimisation();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();
  void FillEventCutFlow(TString cut, TString label , float weight);


  void RunAnalysis(TString plottag, TString tightelid, TString vetoelid, TString looseelid);

  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);  
  float IsDiLep(std::vector<snu::KMuon> muons);

  bool CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets,  TString name, float w);
  bool CheckSignalRegionNN( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KJet> jets, TString name, float w);

  void FillTriggerEfficiency(TString cut, float w, TString label,  std::vector<TString> list);

  float GetTightWeight();
  float  GetMediumWeight();
  void GetSSSignalEfficiency(float w);
  void GetOSSignalEfficiency(float w);
  void GetTriggEfficiency();
  void SignalValidation();
  void RunAnalysis();
  void OptimiseID(bool isss,bool dilep, bool removed0,float w);
  

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
 

  ClassDef ( HNDiMuonOptimisation, 1);
};
#endif
