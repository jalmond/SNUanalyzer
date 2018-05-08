#ifndef GetFakeSF_h
#define GetFakeSF_h

#include "AnalyzerCore.h"

class GetFakeSF : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  GetFakeSF();
  ~GetFakeSF();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();



 private:

 

  ClassDef ( GetFakeSF, 1);
};
#endif
