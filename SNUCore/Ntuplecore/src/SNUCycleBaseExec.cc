#include "SNUCycleBaseExec.h"

ClassImp( SNUCycleBaseExec );

SNUCycleBaseExec::SNUCycleBaseExec() {

}

void SNUCycleBaseExec::BeginCycle()throw( SNUError ){
}

void SNUCycleBaseExec::BeginEvent()throw( SNUError ){

}

void SNUCycleBaseExec::SetUpEvent(Long64_t jentry , float weight)throw( SNUError ){
  m_logger << INFO << "Default SetUpEvent function with "  << jentry << " and weight "<< weight  << SNULogger::endmsg;
}

void SNUCycleBaseExec::ExecuteEvents()throw( SNUError ){
  
}


void SNUCycleBaseExec::EndEvent()throw( SNUError ){

}

void SNUCycleBaseExec::EndCycle()throw( SNUError ){

}


void SNUCycleBaseExec::ClearOutputVectors() throw( SNUError ){

}

void SNUCycleBaseExec::WriteHistograms() throw( SNUError ){

}


