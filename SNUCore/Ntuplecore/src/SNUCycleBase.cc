#include "SNUCycleBase.h"

ClassImp( SNUCycleBase );

SNUCycleBase::SNUCycleBase(){

  m_logger << VERBOSE << "SCycleBase constructed" << SNULogger::endmsg;

}

SNUCycleBase::~SNUCycleBase(){
  m_logger << VERBOSE << "SCycleBase destructed" << SNULogger::endmsg;
}
