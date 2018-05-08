#ifndef EventBase_h
#define EventBase_h

#include <iostream>
using namespace std;


#include "TLorentzVector.h"
#include <vector>
#include "SNUEvent.h"

/// Selection codes
#include "MuonSelection.h"
#include "ElectronSelection.h"
#include "PhotonSelection.h"
#include "JetSelection.h"
#include "FatJetSelection.h"
#include "GenJetSelection.h"
#include "GenSelection.h"
#include "EventSelection.h"

/// SNU codes
#include "KMuon.h"
#include "KElectron.h"
#include "KPhoton.h"
#include "KGenJet.h"
#include "KJet.h"
#include "KFatJet.h"
#include "KTruth.h"
#include "KTrigger.h"

class EventBase {

 public:
  EventBase();
  EventBase(SNUEvent kbase);
  EventBase(EventBase& evbase);
  EventBase& operator= (const EventBase& obj);
  ~EventBase();
  
  inline MuonSelection* GetMuonSel() const {return k_muonsel;}
  inline ElectronSelection* GetElectronSel() const {return k_electronsel;}
  inline PhotonSelection* GetPhotonSel() const {return k_photonsel;}
  inline JetSelection* GetJetSel() const {return k_jetsel;}
  inline FatJetSelection* GetFatJetSel() const {return k_fatjetsel;}
  inline GenJetSelection* GetGenJetSel() const {return k_genjetsel;}
  inline GenSelection* GetTruthSel() const {return k_truthsel;}
  inline EventSelection* GetEventSel() const {return k_eventsel;}

  inline std::vector<snu::KMuon> GetMuons() const {return k_SNUevent->GetMuons();}
  inline std::vector<snu::KElectron> GetElectrons() const {return k_SNUevent->GetElectrons();}
  inline std::vector<snu::KPhoton> GetPhotons() const {return k_SNUevent->GetPhotons();}
  inline std::vector<snu::KJet> GetJets() const {return k_SNUevent->GetJets();}
  inline std::vector<snu::KFatJet> GetFatJets() const {return k_SNUevent->GetFatJets();}
  inline std::vector<snu::KGenJet> GetGenJets() const {return k_SNUevent->GetGenJets();}
  inline std::vector<snu::KTruth> GetTruth() const {return k_SNUevent->GetTruth();}
  inline snu::KTrigger GetTrigger() const {return k_SNUevent->GetTrigger();}
  inline snu::KEvent GetEvent() const {return k_SNUevent->GetEvent();}

  void SetEventBase(snu::KEvent); 

  
  ///Copy constructor
  EventBase(const EventBase& sb);
  
  SNUEvent* GetEventBase() const {return k_SNUevent;}

  SNUEvent* k_SNUevent;
  MuonSelection* k_muonsel;
  ElectronSelection* k_electronsel;
  PhotonSelection* k_photonsel;
  JetSelection* k_jetsel;
  FatJetSelection* k_fatjetsel;
  GenJetSelection* k_genjetsel;
  GenSelection* k_truthsel;
  EventSelection* k_eventsel;
  
 
};

#endif
