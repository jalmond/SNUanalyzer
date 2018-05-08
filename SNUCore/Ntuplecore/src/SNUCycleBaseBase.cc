#include "SNUCycleBaseBase.h"

ClassImp( SNUCycleBaseBase );


SNUCycleBaseBase::SNUCycleBaseBase() : m_logger( "NameNotSet" ) {
  REPORT_VERBOSE( "SNUCycleBaseBase constructed" );
}


void SNUCycleBaseBase::SetLogName( const char* name ) {

  m_logger.SetSource( name );
  return;

}
