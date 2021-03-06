#ifndef SKTreeMakerDiLep_h
#define SKTreeMakerDiLep_h

#include "AnalyzerCore.h"
#include "KEvent.h"
#include "KTrigger.h"


class SKTreeMakerDiLep : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  SKTreeMakerDiLep();
  ~SKTreeMakerDiLep();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );

  void FillCutFlow(TString cut, float weight);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
  std::vector<snu::KPhoton> out_photons;
  std::vector<snu::KJet> out_jets;
  std::vector<snu::KFatJet> out_fatjets;
  std::vector<snu::KGenJet> out_genjets;
  std::vector<snu::KTruth> out_truth;
  snu::KEvent out_event;
  snu::KTrigger out_trigger;

  int nevents;
  int pass_eventcut;
  int pass_vertexcut;

  ClassDef ( SKTreeMakerDiLep, 1);
};
#endif
