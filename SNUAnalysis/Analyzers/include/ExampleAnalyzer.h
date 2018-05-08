#ifndef ExampleAnalyzer_h
#define ExampleAnalyzer_h

#include "AnalyzerCore.h"
class ExampleAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ExampleAnalyzer();
  ~ExampleAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ExampleAnalyzer, 1);
};
#endif
