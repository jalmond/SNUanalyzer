#ifndef SNUCycleBaseBase_H
#define SNUCycleBaseBase_H

// ROOT include(s):
#include <TObject.h>

// Local include(s):
#include "SNULogger.h"

class SNUCycleBaseBase {
 public:
  SNUCycleBaseBase();

  virtual ~SNUCycleBaseBase(){}


 protected:
  /// Function used to set the name of the current cycle
  void SetLogName( const char* name );
  /// Function for accessing the logger object
  SNULogger& logger() const { return m_logger; }

  mutable SNULogger m_logger;

  ClassDef( SNUCycleBaseBase, 0 );

};
#endif
