#include "GenSelection.h"

GenSelection::GenSelection(){
}

GenSelection::GenSelection(SNUEvent ev): BaseSelection(){
  k_snuevent = ev;  
};

GenSelection::~GenSelection() {};

void GenSelection::Selection(std::vector<snu::KTruth>& truthColl) {

  std::vector<snu::KTruth> truth_particles = k_snuevent.GetTruth();  
  
  for (std::vector<snu::KTruth>::iterator trit = truth_particles.begin(); trit!=truth_particles.end(); trit++){
    
    truthColl.push_back(*trit);
  }

  return;
}




GenSelection& GenSelection::operator= (const GenSelection& ms) {
  if(this != &ms){ 
    BaseSelection::operator = (ms);
    k_snuevent = ms.k_snuevent;  
  }
  return *this;
};

GenSelection::GenSelection(const GenSelection& ms):
  BaseSelection(ms)
{
  k_snuevent = ms.k_snuevent;  
};
