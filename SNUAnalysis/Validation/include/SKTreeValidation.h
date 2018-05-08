#ifndef SKTreeValidation_h
#define SKTreeValidation_h

#include "AnalyzerCore.h"

class SKTreeValidation : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  SKTreeValidation();
  ~SKTreeValidation();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();
  void MakeDiMuonValidationPlots(TString muid, float w, float pu_weight,  std::vector<TString> trignames,TString elid, TString jetid, TString tag, bool smearjets, int mcperiod);
  void MakeMuonValidationPlots(TString muid, float w, float pu_weight,  std::vector<TString> trignames,TString elid, TString jetid, TString tag);
  void MakeElectronValidationPlots(TString elid, float w, float pu_weight,  std::vector<TString> trignames,TString muid, TString jetid, TString tag);
  
  void MakeDiElectronValidationPlots(TString elid, float w, float pu_weight,  std::vector<TString> trignames,TString muid, TString jetid, TString tag);
  
  void MakeElMuonValidationPlots(TString id, float w, float pu_weight,  std::vector<TString> trignames,TString elid, TString jetid, TString tag);
  
  void FillCutFlow(TString cut, float weight);
  void counter(TString cut, float w);
  
 private:
  
  std::map<TString, float> mapcounter;
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( SKTreeValidation, 1);
};
#endif
