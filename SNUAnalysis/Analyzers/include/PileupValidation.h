#ifndef PileupValidation_h
#define PileupValidation_h

#include "AnalyzerCore.h"

class PileupValidation : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  PileupValidation();
  ~PileupValidation();

  enum FUNC {
    VALIDATION=0,
    ANALYSIS=1, 
    OPT=2
  };

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );

  void MakeHistograms();


  void counter(TString cut, float w);

 private:

  FUNC functionality ;

  std::map<TString, float> mapcounter;
  
 

  ClassDef ( PileupValidation, 1);
};
#endif
