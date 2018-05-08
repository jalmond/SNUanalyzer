#ifndef GetSignalEff_h
#define GetSignalEff_h

#include "AnalyzerCore.h"

class GetSignalEff : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  GetSignalEff();
  ~GetSignalEff();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();
  void FillEventCutFlow(int cf,TString cut,  float weight);


  bool CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> el,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString name, float w);
  void MatchedJets(std::vector<snu::KJet> jets, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<int> ijets ,TString label);
				      
  
 private:

 

  ClassDef ( GetSignalEff, 1);
};
#endif
