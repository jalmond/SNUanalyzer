#include "SKTreeFiller.h"
#include <stdio.h>  

#include <stdlib.h>
#include <iostream>

using namespace snu;
using namespace std;


SKTreeFiller::SKTreeFiller() :Data() {
  
  TString fitParametersFile = "";
};


SKTreeFiller::~SKTreeFiller() {};


bool SKTreeFiller::SkipTrigger(TString tname){
  
  m_logger << DEBUG << "Trigger: " << tname << SNULogger::endmsg;  
  /// Remove extra unnecisary  triggers (from v7-6-4+ this will not be needed))
  if((tname.Contains("Jpsi")
       || tname.Contains("NoFilters")
       || tname.Contains("Upsilon")
       || tname.Contains("7p5")
       || tname.Contains("Save")
       || tname.Contains("R9Id")
       || tname.Contains("PFMET")
       || tname.Contains("PFHT")
       || tname.Contains("NoHE")
       || tname.Contains("HE10")
       || tname.Contains("PFJet50")
       || tname.Contains("Boost")
       || tname.Contains("LooseIso")
       || tname.Contains("MediumIso")
       || tname.Contains("Mass")
       || tname.Contains("Central")
       || tname.Contains("MW")
       || tname.Contains("EBOnly_VBF")
       || tname.Contains("dEta18"))) return true;
  
  return false;
}


snu::KTrigger SKTreeFiller::GetTriggerInfo(std::vector<TString> trignames){
  snu::KTrigger ktrigger;

  if(!SNUinput){
    ktrigger = *k_inputtrigger;
    return ktrigger;
  }
  m_logger << DEBUG << "Filling trigger Info" << SNULogger::endmsg;


  std::vector<std::string> vHLTInsideDatasetTriggerNames;
  std::vector<bool> vHLTInsideDatasetTriggerDecisions;
  std::vector<int> vHLTInsideDatasetTriggerPrescales;
  

  /// trignames should only be empty id user is running on Catuples and not SKTreeMaker. In this case all triggers are used 
  if(trignames.size() == 0 ){
    for (UInt_t i=0; i< HLT_TriggerName->size(); i++) {
      std::string tgname = HLT_TriggerName->at(i);
      Int_t ps = HLT_TriggerPrescale->at(i);
      vHLTInsideDatasetTriggerNames.push_back(tgname);
      if(ps > 0) vHLTInsideDatasetTriggerDecisions.push_back(true);
      else vHLTInsideDatasetTriggerDecisions.push_back(false);
      vHLTInsideDatasetTriggerPrescales.push_back(ps);
    }
  }

  
  /// vtrigname is vector of ALL triggers in Catuples
  for (UInt_t i=0 ; i< HLT_TriggerName->size(); i++) {
    // trignames is vector of trigger names that we want to store in SKTrees
    // trigname contains names substrings X (where X is for example "HLT_mu") and we store all triggers that start with X

    
    std::string tgname = HLT_TriggerName->at(i);
    if(TString(CatVersion).Contains("v7-6-2")) {
      if(SkipTrigger(TString(tgname)))continue;
    }

    Int_t ps = HLT_TriggerPrescale->at(i);

    for (std::vector<TString>::reverse_iterator it (trignames.end());
	 it != std::vector<TString>::reverse_iterator (trignames.begin());
	 ++it) {

      TString tmpHLT = HLT_TriggerName->at(i);
      if ( tmpHLT.BeginsWith(*it)){
	
	vHLTInsideDatasetTriggerNames.push_back(tgname);
	if(ps > 0) vHLTInsideDatasetTriggerDecisions.push_back(true);
	else vHLTInsideDatasetTriggerDecisions.push_back(false);
	vHLTInsideDatasetTriggerPrescales.push_back(ps);
	
	// if trigger is accepted break from loop
	break;
      }
    } // end of trignames loop
  }// loop of all triggers  
  
  ktrigger.SetHLTInsideDatasetTriggerNames(vHLTInsideDatasetTriggerNames);
  ktrigger.SetHLTInsideDatasetTriggerDecisions(vHLTInsideDatasetTriggerDecisions);
  ktrigger.SetHLTInsideDatasetTriggerPrescales(vHLTInsideDatasetTriggerPrescales);
    
  return ktrigger;
  
}

snu::KEvent SKTreeFiller::GetEventInfo(){
 
  snu::KEvent kevent;

  if(!SNUinput){
    kevent = *k_inputevent;
    if(k_cat_version < 3){
      if(!TString(kevent.CatVersion()).Contains("v7-4"))kevent.SetCatVersion(CatVersion);
    }
    return kevent;
  }
  //  lumimask = snu::KEvent::gold

  m_logger << DEBUG << "Filling Event Info" << SNULogger::endmsg;
  
  // New variable to set catversion. Add this to flat ntuples for next iteration
  kevent.SetCatVersion(CatVersion);

  /// type 1
  // type 1 + ohi corrections
  double met_type1xy = sqrt(pfMET_Type1_PhiCor_Px*pfMET_Type1_PhiCor_Px + pfMET_Type1_Py*pfMET_Type1_Py);
  double phi_type1xy =  TMath::ATan2(pfMET_Type1_Py,pfMET_Type1_Px);
  
  
  /// Default MET is now xy shifted typ1
  if(IsData)  {
    kevent.SetMET(snu::KEvent::pfmet, met_type1xy, phi_type1xy, pfMET_Type1_PhiCor_SumEt);
    kevent.SetPFMETx(pfMET_Type1_PhiCor_Px);
    kevent.SetPFMETy(pfMET_Type1_PhiCor_Py);

    /// Also for completness store type1 without phi corrections
    kevent.SetPFMETType1x(pfMET_Type1_Px);
    kevent.SetPFMETType1y(pfMET_Type1_Py);
    kevent.SetPFMETType1SumEt(pfMET_Type1_SumEt);
    
  }
    /// set unsmeared met variables
  kevent.SetPFMETType1Unsmearedx(pfMET_Type1_Px);
  kevent.SetPFMETType1Unsmearedy(pfMET_Type1_Py);
  kevent.SetPFMETType1xyUnsmearedx(pfMET_Type1_PhiCor_Px);
  kevent.SetPFMETType1xyUnsmearedy(pfMET_Type1_PhiCor_Py);
    
  double topreweight=1.;
  bool settopweight=false;
  if(k_sample_name.Contains("TTLL_powheg"))settopweight=true;
  if(k_sample_name.Contains("TTLJ_powheg"))settopweight=true;
  if(k_sample_name.Contains("TT_powheg"))settopweight=true;
  if(k_sample_name.Contains("TTJets_aMC"))settopweight=true;
  
  
  if(settopweight){
    for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
      if(fabs(gen_PID->at(itx))==6 && fabs(gen_status->at(itx))<30 && fabs(gen_status->at(itx))>20){
	topreweight*=exp(0.0615-0.0005*gen_pt->at(itx));
      }
    }
  }

  kevent.SetTopPtReweight(topreweight);

  if(jet_rho){
    if(jet_rho->size() > 0)kevent.SetRho(jet_rho->at(0));
    else kevent.SetRho(-999.);
  }
  m_logger << DEBUG << "Filling Event Info [2]" << SNULogger::endmsg;
  /// Since some versions of catuples have no metNoHF due to bug in met code 
  
  if(PDFWeights_Error){
    if(PDFWeights_Error->size() > 0){
      std::vector<double>* w1= PDFWeights_Error;
      std::vector<double> w1store;

      for(unsigned int i=0; i < w1->size(); i++){
	w1store.push_back(w1->at(i));
      }
      kevent.SetPDFWeights(w1store);
    }
  }
  if(PDFWeights_Scale){
    if(PDFWeights_Scale->size() > 0){
      std::vector<double>* w1= PDFWeights_Scale;
      std::vector<double> w1store;

      for(unsigned int i=0; i < w1->size(); i++){
        w1store.push_back(w1->at(i));
      }

      kevent.SetScaleWeights(w1store);
    }
  }
  
  
  if(!IsData){
    float jpx(0.), jpy(0.), sjpx(0.), sjpy(0.), sjpxup(0.), sjpxdown(0.),sjpyup(0.), sjpydown(0.) ;
    
    /// only smear jets not close to leptons (use top projection id)
    for(unsigned int ij = 0 ; ij < jet_pt->size(); ij++){
      bool close_to_lepton(false);
      if(jet_pt->at(ij) < 10.) continue;
      for(unsigned int im=0; im < muon_pt->size(); im++){
	if(muon_pt->at(im) < 10.) continue;
	if(fabs(muon_eta->at(im)) > 2.5) continue;
	// find full definition for 13 TeV
	//if(muon_relIso04->at(im) > 0.2)  continue;
        double dr = sqrt( pow(fabs( jet_eta->at(ij) - muon_eta->at(im)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( jet_phi->at(ij) - muon_phi->at(im))),2.0));
	if(dr < 0.4){
	  close_to_lepton=true;
	}
      }
      for(unsigned int iel=0; iel < electron_pt->size(); iel++){
	if(electron_pt->at(iel) < 10.) continue;
        if(fabs(electron_eta->at(iel)) > 2.5) continue;
	// find full definition for 13 TeV                                                                                                                                          if(electrons_relIso03->at(ilep) > 0.15)  continue;
        double dr = sqrt( pow(fabs( jet_eta->at(ij) - electron_eta->at(iel)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( jet_phi->at(ij) - electron_phi->at(iel))),2.0));
        if(dr < 0.4){
          close_to_lepton=true;
        }
      }
      
      if(close_to_lepton) continue;
      
      float jet_px = jet_pt->at(ij) *TMath::Cos(jet_phi->at(ij)); 
      float jet_py = jet_pt->at(ij) *TMath::Sin(jet_phi->at(ij));
      jpx +=  jet_px;
      jpy +=  jet_py;
      
      // JS jet_smearedRes ?
      if(!IsData){
	sjpx +=  1. *jet_px;
	sjpy +=  1. *jet_py;
      }
      else{
	sjpx +=  jet_px;
        sjpy +=  jet_py;
      }
      sjpxup +=  1. *jet_px;
      sjpyup +=  1. *jet_py;
      
      sjpxdown +=  1. *jet_px;
      sjpydown +=  1. *jet_py;

    }

    // met_jetRes_Px_up ==met_Px since no smearing is applied in miniaods -> cattools
    float met_x  = pfMET_Type1_PhiCor_Px  +  jpx - sjpx;
    float met_y  = pfMET_Type1_PhiCor_Py  +  jpy - sjpy;
    float met_newpt = sqrt(met_x*met_x+ met_y*met_y);
    float met_newphi = TMath::ATan2(met_y,met_x);
    
    kevent.SetMET(snu::KEvent::pfmet,  met_newpt,met_newphi, pfMET_Type1_PhiCor_SumEt);  
    kevent.SetPFMETx(met_x);
    kevent.SetPFMETy(met_y);

    /// correct MET for jets smearing
    float type1_met_x  = pfMET_Type1_Px  +  jpx - sjpx;
    float type1_met_y  = pfMET_Type1_Py +   jpy - sjpy;
    
    kevent.SetPFMETType1x(type1_met_x);
    kevent.SetPFMETType1y(type1_met_y);		  
    kevent.SetPFMETType1SumEt(pfMET_Type1_SumEt);           

    /// Fix met phi 
    float met_x_jer_up  = pfMET_Type1_PhiCor_Px +  jpx - sjpxup;
    float met_y_jer_up   = pfMET_Type1_PhiCor_Py  +  jpy - sjpyup;
    float met_newpt_jerup = sqrt(met_x_jer_up*met_x_jer_up+ met_y_jer_up*met_y_jer_up);
    float met_x_jer_down   = pfMET_Type1_PhiCor_Px  +  jpx - sjpxdown;
    float met_y_jer_down  = pfMET_Type1_PhiCor_Py  +  jpy -sjpydown;
    float met_newpt_jerdown = sqrt(met_x_jer_down*met_x_jer_down+ met_y_jer_down*met_y_jer_down);

      
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::JetRes,     met_newpt_jerup);
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::JetRes,     met_newpt_jerdown);
    
  }
  

  m_logger << DEBUG << "Filling Event Info [3]" << SNULogger::endmsg;
  

  /*
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::MuonEn,     sqrt(met_muonEn_Px_up*met_muonEn_Px_up + met_muonEn_Py_up*met_muonEn_Py_up));
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::MuonEn,     sqrt(met_muonEn_Px_down*met_muonEn_Px_down + met_muonEn_Py_down*met_muonEn_Py_up));
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::ElectronEn, sqrt(met_electronEn_Px_up*met_electronEn_Px_up + met_electronEn_Py_up*met_electronEn_Py_up));
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::ElectronEn, sqrt(met_electronEn_Px_down*met_electronEn_Px_down + met_electronEn_Py_down*met_electronEn_Py_down));
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::Unclustered,sqrt(met_unclusteredEn_Px_up->at(0)*met_unclusteredEn_Px_up->at(0) + met_unclusteredEn_Py_up->at(0)*met_unclusteredEn_Py_up->at(0)));
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::Unclustered,sqrt(met_unclusteredEn_Px_down->at(0)*met_unclusteredEn_Px_down->at(0) + met_unclusteredEn_Py_down->at(0)*met_unclusteredEn_Py_down->at(0)));
    kevent.SetPFSumETShift(snu::KEvent::up,     snu::KEvent::Unclustered,met_unclusteredEn_SumEt_up->at(0));
    kevent.SetPFSumETShift(snu::KEvent::down,   snu::KEvent::Unclustered,met_unclusteredEn_SumEt_down->at(0));
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::JetEn,      sqrt(met_jetEn_Px_up->at(0)*met_jetEn_Px_up->at(0) + met_jetEn_Py_up->at(0)*met_jetEn_Py_up->at(0)));
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::JetEn,      sqrt(met_jetEn_Px_down->at(0)*met_jetEn_Px_down->at(0) + met_jetEn_Py_down->at(0)*met_jetEn_Py_down->at(0)));
    kevent.SetPFSumETShift(snu::KEvent::up,     snu::KEvent::JetEn,      met_jetEn_SumEt_up->at(0));
    kevent.SetPFSumETShift(snu::KEvent::down,   snu::KEvent::JetEn,      met_jetEn_SumEt_down->at(0));
  */
  /// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/python/patPFMETCorrections_cff.py
  /// jets > 15 GeV in mc smeared. This is not done in cattolls so branches have no change,
  /// Apply this here
  /// 
  
  
  
  m_logger << DEBUG << "Filling Event Info [4]" << SNULogger::endmsg;
  
  /// Filling event variables
    
  kevent.SetIsData(IsData);
  kevent.SetRunNumber(run);
  kevent.SetEventNumber(event);
  kevent.SetLumiSection(lumi);
  
  if(!IsData){
  
    kevent.SetPUWeight(snu::KEvent::central,double(PUweight)); // JSKIM
    kevent.SetPUWeight(snu::KEvent::down,double(pileUpReweightMinus));
    kevent.SetPUWeight(snu::KEvent::up,  double(pileUpReweightPlus));

  }
  kevent.SetGenId(genWeight_id1, genWeight_id2);
  // kevent.SetLHEWeight(lheWeight);   JSKIM
  kevent.SetGenX(genWeight_X1, genWeight_X2);
  kevent.SetGenQ(genWeight_Q);
  if(gen_weight > 0.) kevent.SetWeight(1.);
  else kevent.SetWeight(-1.);
  
  kevent.SetVertexInfo(vertex_X, vertex_Y, vertex_Z,0. );
  
  /// MET filter cuts/checks

  
  /// 
  //  kevent.SetPileUpInteractionsTrue(nTrueInteraction); // JSKIM
    
  kevent.SetNVertices(nPV);
  kevent.SetNGoodVertices(nPV); // JSKIM
  
  //  kevent.SetIsGoodEvent(nGoodPV);

  /// MET filter cuts/checks
  kevent.SetPassEcalDeadCellTriggerPrimitiveFilter(Flag_EcalDeadCellTriggerPrimitiveFilter);
  kevent.SetPassHBHENoiseFilter(Flag_HBHENoiseFilter);
  kevent.SetPassHBHENoiseIsoFilter(Flag_HBHENoiseIsoFilter);

  kevent.SetPassTightHalo2016Filter(Flag_globalTightHalo2016Filter);
  kevent.SetPassBadChargedCandidateFilter(Flag_BadChargedCandidateFilter);
  kevent.SetPassBadPFMuonFilter(Flag_BadPFMuonFilter);
  
  // JSKIM

  return kevent;
}


std::vector<KPhoton> SKTreeFiller::GetAllPhotons(){

  std::vector<KPhoton> photons;

  
  if(!SNUinput){
    for(std::vector<KPhoton>::iterator kit  = k_inputphotons->begin(); kit != k_inputphotons->end(); kit++){
      photons.push_back(*kit);
    }
    return photons;
  }
  for (UInt_t iph=0; iph< photon_eta->size(); iph++) {
    if(photon_pt->at(iph) != photon_pt->at(iph)) continue;
    KPhoton ph;
    
    ph.SetPtEtaPhiM(photon_pt->at(iph),photon_eta->at(iph), photon_phi->at(iph),0.); // JSKIM 

    ph.SetIsLoose(photon_passLooseID->at(iph));
    ph.SetIsMedium(photon_passMediumID->at(iph));
    ph.SetIsTight(photon_passTightID->at(iph));
    ph.SetPassMVA(photon_passMVAID_WP80->at(iph));
    //ph.SetMCMatched(photon_mcMatched->at(iph)); 
    //ph.SetHasPixSeed(photon_haspixseed->at(iph));
    //ph.SetPassElVeto(photon_passelectronveto->at(iph));

    ph.SetChargedHadIsoNoEA(photon_ChIso->at(iph));
    //ph.SetpuChargedHadIsoNoEA();  JSKIM
    ph.SetNeutalHadIsoNoEA(photon_NhIso->at(iph));
    ph.SetPhotonIsoNoEA(photon_PhIso->at(iph));
    //ph.SetRhoIso(); JSKIM
    ph.SetChargedHadIso(photon_ChIsoWithEA->at(iph));
    ph.SetPhotonIso(photon_PhIsoWithEA->at(iph));
    ph.SetNeutalHadIso(photon_NhIsoWithEA->at(iph));
    ph.SetSigmaIetaIeta(photon_Full5x5_SigmaIEtaIEta->at(iph));
    //ph.SetR9(->at(iph)); JSKIM
    ph.SetHoverE(photon_HoverE->at(iph));
    ph.SetSCEta(photon_scEta->at(iph));
    ph.SetSCPhi(photon_scPhi->at(iph));
    //ph.SetSCRawE(->at(iph)); JSKIM
    //ph.SetSCPreShowerE(->at(iph)); JSKIM
    
    photons.push_back(ph);
  }
  std::sort( photons.begin(), photons.end(), isHigherPt );

  return photons;

}

std::vector<KElectron> SKTreeFiller::GetAllElectrons(){

  std::vector<KElectron> electrons;

  if(!SNUinput){
    for(std::vector<KElectron>::iterator kit  = k_inputelectrons->begin(); kit != k_inputelectrons->end(); kit++){
      electrons.push_back(*kit);
    }
    return electrons;
  }

  m_logger << DEBUG << "Filling electron Info " << electron_eta->size() << SNULogger::endmsg;
  
  vector<int> matched_truth;
  for (UInt_t iel=0; iel< electron_eta->size(); iel++) {
    
    if(electron_pt->at(iel) != electron_pt->at(iel))    continue;
    
    KElectron el;

    /// Kinematic Variables
    el.SetPtEtaPhiE(electron_pt->at(iel),electron_eta->at(iel), electron_phi->at(iel),electron_Energy->at(iel));

    el.SetSmearFactor(electron_Energy_Smear_Up->at(iel)/ electron_Energy->at(iel));
    //el.SetTrigMatch(electron_trigmatch->at(iel)); JSKIM
    el.SetSCEta(electron_scEta->at(iel));
   
    el.Setdz( electron_dz->at(iel));
    el.Setdxy(electron_dxy->at(iel) );
    if(electron_mva){
      el.SetMVA(electron_mva->at(iel) );
      el.SetZZMVA(electron_zzmva->at(iel) );
    }
    if(electron_sigdxy){
      if(electron_sigdxy->size() > 0 )el.Setdxy_sig(electron_sigdxy->at(iel) );
    }
    el.SetPFChargedHadronIso(0.3, electron_puChIso03->at(iel));
    el.SetPFPhotonIso(0.3,electron_phIso03->at(iel));
    el.SetPFNeutralHadronIso(0.3,electron_nhIso03->at(iel));
    el.SetPFRelIsoRho(0.3,electron_relIsoRho03->at(iel));
    el.SetPFRelIsoBeta(0.3,electron_relIsoBeta03->at(iel));



    m_logger << DEBUG << "Filling electron_minirelIso " << SNULogger::endmsg;
    if(electron_minirelIso) el.SetPFRelMiniIso(electron_minirelIso->at(iel));
    
    m_logger << DEBUG << "Filling electron Info 2" << SNULogger::endmsg;
    
    el.SetPFChargedHadronIso(0.4,electron_puChIso04->at(iel));
    el.SetPFPhotonIso(0.4,electron_phIso04->at(iel));
    el.SetPFNeutralHadronIso(0.4,electron_nhIso04->at(iel));
    el.SetPFRelIso(0.4,electron_relIso04->at(iel));
    
    el.SetPFAbsIso(0.3,electron_absIso03->at(iel));
    el.SetPFAbsIso(0.4,electron_absIso04->at(iel));


    /// set Charge variables
    el.SetCharge(electron_q->at(iel));
    el.SetGsfCtfScPixCharge(electron_isGsfCtfScPixChargeConsistent->at(iel));
    
    m_logger << DEBUG << "Filling electron Info 3" << SNULogger::endmsg;
    /// set conversion variables
    
    if(electron_shiftedEnDown){
      el.SetShiftedEUp(electron_shiftedEnUp->at(iel));
      el.SetShiftedEDown(electron_shiftedEnDown->at(iel));
    }

    if(electron_missinghits)el.SetMissingHits(electron_missinghits->at(iel));
    el.SetSNUID(electron_electronID_snu->at(iel));
    el.SetPassVeto(electron_electronID_veto->at(iel));
    el.SetPassLoose(electron_electronID_loose->at(iel));
    el.SetPassMedium(electron_electronID_medium->at(iel));
    el.SetPassTight(electron_electronID_tight->at(iel));
    if(electron_electronID_hlt)el.SetPassHLT(electron_electronID_hlt->at(iel));
    /// HEEP
    //el.SetPassHEEP(electron_electronID_heep->at(iel));

    // MVA
    el.SetPassMVATrigMedium(electron_electronID_mva_trig_medium->at(iel));
    el.SetPassMVATrigTight(electron_electronID_mva_trig_tight->at(iel));
    el.SetPassMVANoTrigMedium(electron_electronID_mva_medium->at(iel));
    el.SetPassMVANoTrigTight(electron_electronID_mva_tight->at(iel));
    if(electron_electronID_mva_zz)el.SetPassMVAZZ(electron_electronID_mva_zz->at(iel));

    el.SetIsPF(electron_isPF->at(iel));
    if(electron_isTrigMVAValid) el.SetIsTrigMVAValid(electron_isTrigMVAValid->at(iel));
    //el.SetIsMCMatched(electron_mcMatched->at(iel));
    el.SetHasMatchedConvPhot(electron_passConversionVeto->at(iel));
    
    el.SetTrkVx(electron_x->at(iel));
    el.SetTrkVy(electron_y->at(iel));
    el.SetTrkVz(electron_z->at(iel));
    m_logger << DEBUG << "Filling electron Info 4" << SNULogger::endmsg;    

    //// Set Is ChargeFlip
    bool isprompt= false;
    bool from_tau = false;
    
    int mother_index=-1;
    int mother_pdgid=-1;
    int matched_index=-1;
    int mc_pdgid=-1;
    bool matched_in_Dr=false;

    int           eltype=0;
    bool conv_veto=false;
    if(k_cat_version  > 3){
      
      if(gen_pt){
	// Default deltaR setting for matching
	float min_Dr=0.1;
	/// Loop over all gen particles
	for (UInt_t it=0; it< gen_pt->size(); it++ ){
	  
	  
	  /// Requirements to make sure no crash or warnings with pt=0
	if(gen_mother_index->at(it) <= 0)continue;
	if(gen_mother_index->at(it) >= int(gen_pt->size()))continue;
	if(gen_pt->at(it) < 0.001) continue;
	

	double match_eta =electron_eta->at(iel);
	double match_pt =electron_pt->at(iel);
	double match_phi =electron_phi->at(iel);
	double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));

	
	
	/// check for photon conversion veto in DY/ZG                                                                                                                               
	if(fabs(gen_PID->at(it)) ==22){
	  if(gen_pt->at(it) > 10.){	
	    if(dr < 0.3){
	      if(gen_isprompt->at(it) && gen_status->at(it) ==1) {
		conv_veto=true;
		for (UInt_t it_ph=0; it_ph< gen_pt->size(); it_ph++ ){
		  if(it==it_ph) continue;

		  // check ph is matched to q or g from matrix element (st 23)
		  if(gen_status->at(it_ph) ==23){
		    if(fabs(gen_PID->at(it_ph)) < 7 || fabs(gen_PID->at(it_ph))==21){
		      double drph = sqrt( pow(fabs(gen_eta->at(it_ph) - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( gen_phi->at(it_ph) - gen_phi->at(it))),2.0));
		      if(drph < 0.05) conv_veto=false;
		    }
		  }
		}
	      }	      
	    }
	  }
	}

	/// Matching using instructions on
	/// https://indico.cern.ch/event/292928/contributions/1650088/attachments/547844/755123/talk_electron_contribution.pdf
	/// 

	/// Match required to status 1 electron
	if(gen_status->at(it) != 1) continue;
	if(fabs(gen_PID->at(it)) != 11) continue;
	if(gen_pt->at(it) < 5.) continue;

	/// Check status 1 electron is not matched already to areco electron
	bool already_matched=false;
	for(unsigned int im=0; im < matched_truth.size();im++){
          if(it == unsigned(matched_truth.at(im))) already_matched=true;
        }
        //if(already_matched) continue;
	  
	if(matched_in_Dr){
	  /// This is for multiple matched status 1 el.
	  /// In case multiple status 1 electrons are matched with same mother check pt
	  if(gen_mother_index->at(it) == gen_mother_index->at(matched_index)){
	    if (dr < 0.1){
	      // in case 2+ electrons from same mother electron (conversion) also match in pt
	      if( fabs(gen_pt->at(it)-electron_pt->at(iel)) < fabs(gen_pt->at(matched_index)-electron_pt->at(iel))) matched_index=it;
	    }
	  }
	  else if ((dr < min_Dr) ){
	    
	    /// find closest match in dR to status 1
            if( (fabs(match_pt - gen_pt->at(it))/gen_pt->at(it)) < 2.) {
	      
	      matched_in_Dr=true;
	      min_Dr= dr;
	      
	      /// set index of matched status 1 electron
	      matched_index=it;
	    }
	  }
	}
	else{
	  /// first match status 1 electron
	  if ((dr < min_Dr) ){
	    if( (fabs(match_pt - gen_pt->at(it))/gen_pt->at(it)) < 2.) {
	      
	      /// find closest match in dR to status 1
	      matched_in_Dr=true;
	      min_Dr= dr;
	      
	      /// set index of matched status 1 electron
	      matched_index=it;
	    }
	  }
	}
      }// end of gen loop to find status 1 electron
	
	//cout << iel << " " << electron_pt->at(iel) << " " << electron_eta->at(iel) << " " << electron_phi->at(iel) << " " << conv_veto << endl;
	
      ///// treat case where there is a matched status 1 electron:
      //// classify into prompt:Fake:FromTau

      if(matched_in_Dr){
	/// Find closest non electron ancesteror
	float pdgid = gen_PID->at(matched_index);

	// mc_pdgid = closest matched status 1 pdgid
	mc_pdgid= int(pdgid);

	// mindex = mother index: will loop to find first non el mother
	int mindex= matched_index;

	while ( (fabs(gen_PID->at(mindex)) == 11)) {
	  pdgid = gen_PID->at(mindex);
	  mindex=gen_mother_index->at(mindex);
	}

	/// pdgid is now of electron from non electron mother
	//  mindex = index for mother of non electron ancestor
	
	if( (fabs(gen_PID->at(mindex)) == 23) || (fabs(gen_PID->at(mindex)) == 24)) {
	  /// Check if el from Z/W is CF and if it is from a photon conversion
	  
	  eltype=1;
	  int n_el_from_el=0;
	  float charge_sum=0.;
	  /// Loop over electrons: Find mother of matched status 1 and see what other daughters there are:
	  /// In case of a conversion i.e  Z->ee->eephoton->eeee the status 23 electorn decays to 3 electrons e+e+e- or e-e-e+
	  bool isthirdel_fromconv(false);
	  for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	    if(gen_mother_index->at(itx) <= 0)continue;
	    if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(itx) < 0.001) continue;
	    if(fabs(gen_PID->at(itx)) ==11) {
	      if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) { 
		charge_sum+= gen_PID->at(itx); n_el_from_el++;
		if(n_el_from_el==3){
		  if(itx == matched_index) isthirdel_fromconv=true;
		}
		if(n_el_from_el==5){
                  if(itx == matched_index) isthirdel_fromconv=true;
		}
	      }

	    }
	  }

	  /// Set if conversion i.ei e->eee
	  /// Two methods: 
	  /// 1) check pdgid of status 1 el vs mother. if < 0 it is a converison
	  /// 2) In case closest status 1 el is not opposite charge truth check number of electrons from mother if 3 it is a conversion
	  if((gen_PID->at(matched_index)  * pdgid) < 0 )  {el.SetIsPhotonConversion(true);           eltype=2;} 
	  else  el.SetIsPhotonConversion(false);
	  
	  
	  if(!isthirdel_fromconv){
	    if(n_el_from_el ==3&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true); }
	    if(n_el_from_el ==5&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true); }
	  }
	  else{
	    if(pdgid * electron_q->at(iel) > 0 )  {
	      if(n_el_from_el ==3&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true);}
	      if(n_el_from_el ==5&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true);}
	    }
	  }

	  /// Check if it is a chargeflip.
	  /// Either from a conversion or just reconstructed charge is wrong
	  if(pdgid * electron_q->at(iel) > 0 )    
	    { el.SetIsChargeFlip(true); 
	      if(eltype == 2 || eltype == 3){
		if(eltype == 2) eltype=4;
		if(eltype == 3) eltype=5;
	      }
	      else eltype=6;
	    }

	  else     el.SetIsChargeFlip(false);
	  
	  mother_index=mindex;
	  mother_pdgid=gen_PID->at(mindex);
	  isprompt=true; /// means is prompt
	  
	}/// end of Z/W
	else {
	  if(gen_status->at(mindex) == 2){
	    if(fabs(gen_PID->at(mindex)) > 50) {isprompt=false; mother_pdgid=gen_PID->at(mindex); mother_index=mindex; from_tau=false;
	      eltype=7;
	      
	      if(gen_isprompt->at(matched_index)){
		cout << "matched FAKE, but isPrompt flag??" << endl;
		cout << "------------------CF "<< endl;
		cout << "gen_isprompt = " << gen_isprompt->at(matched_index)  << endl;
		cout << "gen_isdecayedleptonhadron = " <<gen_isdecayedleptonhadron->at(matched_index)  << endl;
		cout << "gen_isdirecthadrondecayproduct  = " <<gen_isdirecthadrondecayproduct->at(matched_index)  << endl;
		cout << "gen_ishardprocess  = " << gen_ishardprocess->at(matched_index)  << endl;
		cout << "gen_istaudecayproduct =  " << gen_istaudecayproduct->at(matched_index)  << endl;
		cout << "gen_isprompttaudecayproduct =  " <<  gen_isprompttaudecayproduct->at(matched_index)  << endl;
	      }
	      
	    }
	    else {
	      isprompt=true;
	      mother_pdgid=gen_PID->at(mindex); mother_index=mindex; from_tau=false;
	      eltype=8;
	      
	      
	      if(fabs(gen_PID->at(mindex)) == 22){
		if(fabs(gen_PID->at(gen_mother_index->at(mindex))) > 50){
		  eltype=9;
		}
		else {
		  eltype=10;
		}
	      }
	    }
	    
	    if(fabs(gen_PID->at(mindex)) == 15){
	      eltype=11;

	      isprompt=true; mother_pdgid=gen_PID->at(mindex);  mother_index=mindex; from_tau=true;
	      // Check if el from tau  is CF
	      
	      int n_el_from_el=0;
	      float charge_sum=0.;
	      for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		if(itx == matched_index) continue;
		if(gen_mother_index->at(itx) <= 0)continue;
		if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
		if(gen_pt->at(itx) < 0.001) continue;
		if(fabs(gen_PID->at(itx)) ==11) {
		  if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) { charge_sum+= gen_PID->at(itx); n_el_from_el++;
		  }
		}
	      }// end of truth loop to check Conv
	      if((gen_PID->at(matched_index)  * pdgid) < 0 )  el.SetIsPhotonConversion(true);
	      else  el.SetIsPhotonConversion(false);
	      if(n_el_from_el ==3&& (fabs(charge_sum) == 11))  el.SetIsPhotonConversion(true);
	      if(n_el_from_el ==5&& (fabs(charge_sum) == 11))  el.SetIsPhotonConversion(true);

	      if(fabs(gen_PID->at(gen_mother_index->at(mother_index))) > 50) {isprompt=false; eltype=12;}
	      
	      if(pdgid * electron_q->at(iel) > 0 )     {el.SetIsChargeFlip(true); eltype=13; }
	      else     el.SetIsChargeFlip(false);
	    }
	    
	  }/// end of status 2 check
	  else {
	    /// using new method for matching: These events are set as prompt 
	    isprompt=true;mother_pdgid=-99999; mother_index=mindex; from_tau=false; 
	    eltype=14;
	    int n_el_from_eg=0;  
	    vector<KTruth> vel_tmp;
	    bool isthirdel_fromconv(false);
	    bool neutrino_invertex(false);
	    for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	      if(gen_mother_index->at(itx) <= 0)continue;
	      if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
	      if(gen_pt->at(itx) < 0.001) continue;
	      if(fabs(gen_PID->at(itx)) ==11) {
		if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) {  n_el_from_eg++;
		  if(n_el_from_eg==3){isthirdel_fromconv=true; }
		  if(n_el_from_eg==5){isthirdel_fromconv=true; }
		  if(gen_status->at(itx) ==1){
		    KTruth truthe;
		    truthe.SetPtEtaPhiE(gen_pt->at(itx), gen_eta->at(itx), gen_phi->at(itx), gen_energy->at(itx));
		    vel_tmp.push_back(truthe);
		  }
		}
	      }
	      
	      if(fabs(gen_PID->at(itx)) ==12) {

		int index_mother_nu=gen_mother_index->at(itx);
		while (fabs(index_mother_nu) ==12){
		  index_mother_nu=gen_mother_index->at(index_mother_nu);
		}
		
		if(index_mother_nu == mindex) {
		  neutrino_invertex=true;
		}
	      }
	    } // end of truth loop to check Conv
	    
	    if(neutrino_invertex) eltype=15;

	    if(vel_tmp.size() ==2) {
	      KParticle ll = vel_tmp[0] + vel_tmp[1];
	      if(fabs(ll.M()) < 5.) eltype=16;
	    }


	    if((gen_PID->at(matched_index)  * pdgid) < 0 )  {el.SetIsPhotonConversion(true);  eltype=17;}
	    else el.SetIsPhotonConversion(false);
	    
	    if(n_el_from_eg ==3&&!isthirdel_fromconv)  {el.SetIsPhotonConversion(true); eltype=18;}
	    if(isthirdel_fromconv&&n_el_from_eg ==3){
	      if(pdgid * electron_q->at(iel) > 0 )   {
		el.SetIsPhotonConversion(true);
		eltype=18;
	      }
	    }
	    if(isthirdel_fromconv&&n_el_from_eg ==5){
              if(pdgid * electron_q->at(iel) > 0 )   {
                el.SetIsPhotonConversion(true);
                eltype=18;
              }
            }

	    if(pdgid * electron_q->at(iel) > 0 )  {
	      el.SetIsChargeFlip(true);
	      if(eltype==17  || eltype == 18){
		if(eltype==17 ) eltype=19;
		if(eltype==18 ) eltype=20;
	      }
	      else eltype=21;
	    }
	    else     el.SetIsChargeFlip(false);
	    
	    /// speacial treatment for signal 

	    if( fabs(gen_PID->at(mindex))>= 9900012 &&  fabs(gen_PID->at(mindex)) < 9900025 )mother_pdgid= gen_PID->at(mindex);
	    	    
	    
	  }  // not gen status 2
	} // not Z/W daughter
      }  /// In case no status 1 electron is found : classify electron fake
      else{
	if(gen_pt){
	  for (UInt_t it=0; it< gen_pt->size(); it++ ){
	    if(gen_mother_index->at(it) <= 0)continue;
	    if(gen_mother_index->at(it) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(it) < 0.001) continue;
	    
	    double match_eta =electron_eta->at(iel);
	    double match_phi =electron_phi->at(iel);
	    double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));
	  
	    bool already_matched=false;
	    for(unsigned int im=0; im < matched_truth.size();im++){
	      if(it == unsigned(matched_truth.at(im))) already_matched=true;
	    }
	    
	    // in coversion case  2 matched electrons to one photon
	    if(fabs(gen_PID->at(it)) != 22 && already_matched) continue;
	    
	    if (dr <0.1){
	      matched_in_Dr=true;
	      int mindex= gen_mother_index->at(it);
	      float pdgid = gen_PID->at(it);

	      
	      /// Unlikely to have mother as electron but just in case
	      while ( (fabs(gen_PID->at((mindex))) == 11)) {
		mindex=gen_mother_index->at(mindex);
	      }
	      // isprompt = false since it failed status 1 matching
	      isprompt=false;
	      /// mother index of first non electron
	      int grandmother =0;
	      if(gen_mother_index->at(mindex) > 0)grandmother=gen_PID->at(gen_mother_index->at(mindex));
	      mother_pdgid=gen_PID->at(mindex);
	      mother_index=mindex;
	      matched_index = it;
	      mc_pdgid= int(gen_PID->at(it));
	      if(fabs(pdgid) == 22) {

		if(fabs(mother_pdgid) > 50) eltype=22;
		else eltype=23;		
		
		//// This case is not a conversion
		//el.SetIsPhotonConversion(true);
		
		//if(gen_PID->at(gen_mother_index->at(it)) * electron_q->at(iel) > 0 )     el.SetIsChargeFlip(true);
		//else     el.SetIsChargeFlip(false);
		
		from_tau=false;
		break;
	      }
	      
	      else if(fabs(pdgid) == 15){
		from_tau=true;
		if(fabs(mother_pdgid) > 50) eltype=24;
		else eltype=25;
		if(fabs(mother_pdgid)==15){
		  if(fabs(grandmother) > 50) eltype=-24;
		  else eltype=-25;
		}
	      }
	      else if(fabs(pdgid) == 1){
		eltype=26;
	      }
	      else if(fabs(pdgid) == 2){
                eltype=27;
              }
	      else if(fabs(pdgid) == 3){
                eltype=28;
              }
	      else if(fabs(pdgid) == 4){
                eltype=29;
              }
	      else if(fabs(pdgid) == 5){
                eltype=30;
              }
	      else if(fabs(pdgid) == 21){
                eltype=31;
              }
	      else if(fabs(pdgid) == 211){
		eltype=32;
              }
	      else if(fabs(pdgid) == 310){
		eltype=33;
              }
	      else if(fabs(pdgid) == 431){
		eltype=34;
              }
	      else if(fabs(pdgid) == 13){
		eltype=35;
              }
	      else if(fabs(pdgid) == 311){
                eltype=36;
              }
	      else if(fabs(pdgid) > 500 && fabs(pdgid) < 600){
                eltype=37;
              }
	      else eltype=38;
	      
	    }// dr req
	  }// loop over gen vector
	}// require gen info
      }// no status 1 match
      }/// END OF TRUTH MATCHING
      
      matched_truth.push_back(matched_index);
      ///matched_index is index which matches reco muon with smallest dR
    ///- If multiple status 1 muons are matched look at closest in pt
    ///- In no status 1 is matched set as not prompt butlook for closest particle in dR
    /// - In noparticles within dR < 0.1 matched_in_Dr= false
      
      el.SetType(eltype);
    if(!matched_in_Dr){
      el.SetIsMCMatched(false);
      el.SetIsFromTau(false);
      el.SetMCMatchedPdgId(-1);
      el.SetMotherPdgId(0);
      el.SetMotherTruthIndex(-1);
      el.SetMCTruthIndex(-1);
      el.SetIsMCExternalConversion(conv_veto);
      if(conv_veto)el.SetType(40);
    }
    else{
      
      if(!isprompt){
	if((gen_isprompt->at(matched_index) ==1 )&& (gen_status->at(matched_index) == 1)){
	  
	  //cout << "gen_istaudecayproduct =  " << gen_istaudecayproduct->at(matched_index)  << endl;
	  //cout << "gen_isprompttaudecayproduct =  " <<  gen_isprompttaudecayproduct->at(matched_index)  << endl;
	  if(!(gen_istaudecayproduct->at(matched_index)   || gen_isprompttaudecayproduct->at(matched_index))){
	    //cout << "matched as prompt yet status flag is not prompt" << endl;
	    //cout << "matched_index = " << matched_index << endl;
	    //cout << "reco "<< electron_pt->at(iel)<< " " << electron_eta->at(iel)  << " " << electron_phi->at(iel) << endl;;
	    //for (UInt_t it=0; it< gen_pt->size(); it++ ){
	    //	      if(gen_mother_index->at(it) <= 0)continue;
	    //if(gen_mother_index->at(it) >= int(gen_pt->size()))continue;
	    //if(gen_pt->at(it) < 0.1) continue;
	    //cout << it << " " << gen_pt->at(it)  << " " << gen_eta->at(it) << " " << gen_phi->at(it)<< " " << gen_PID->at(it) << "  " << gen_status->at(it) << " " << gen_PID->at(gen_mother_index->at(it)) <<" "  <<  gen_mother_index->at(it) << " " << gen_isprompt->at(it)  <<endl;
	      
	    //}
	  }
	}
      }
      
      //if(gen_isprompt->at(matched_index) 
      el.SetIsMCExternalConversion(conv_veto);
      el.SetIsMCMatched(isprompt);
      el.SetIsFromTau(from_tau);
      el.SetMotherPdgId(mother_pdgid);
      el.SetMCMatchedPdgId(mc_pdgid);
      el.SetMotherTruthIndex(mother_index);
      el.SetMCTruthIndex(matched_index);
      if(gen_status->at(matched_index)==1)el.SetIsPromptFlag(gen_isprompt->at(matched_index));
      if(conv_veto)el.SetType(40);

    }
    }
    electrons.push_back(el);
  }
  m_logger << DEBUG << "END electrons " << SNULogger::endmsg;
  std::sort( electrons.begin(), electrons.end(), isHigherPt );
  
  return electrons;
}


void SKTreeFiller::ERRORMessage(TString comment){
  
  m_logger << ERROR << "SKTreeFiller had a probleming filling " << comment << ". This variable is not present in the current SNUntuples." << SNULogger::endmsg;   
}



std::vector<KGenJet> SKTreeFiller::GetAllGenJets(){

  std::vector<KGenJet> genjets;
  if(IsData) return genjets;
  if(!SNUinput){
    if(k_inputgenjets){
      for(std::vector<KGenJet>::iterator kit  = k_inputgenjets->begin(); kit != k_inputgenjets->end(); kit++){
	genjets.push_back(*kit);
      }
    }
    return genjets;
  }
  if(k_cat_version < 3){
    for (UInt_t ijet=0; ijet< slimmedGenJet_pt->size(); ijet++) {
      KGenJet jet;
      jet.SetPtEtaPhiE(slimmedGenJet_pt->at(ijet), slimmedGenJet_eta->at(ijet), slimmedGenJet_phi->at(ijet), slimmedGenJet_energy->at(ijet));
      genjets.push_back(jet);
    }
    return genjets;
  }
  
  for (UInt_t ijet=0; ijet< genjet_pt->size(); ijet++) {
    KGenJet jet;
    jet.SetPtEtaPhiE(genjet_pt->at(ijet), genjet_eta->at(ijet), genjet_phi->at(ijet), genjet_energy->at(ijet));
    jet.SetGenJetEMF(genjet_emf->at(ijet));
    jet.SetGenJetHADF(genjet_hadf->at(ijet));
    jet.SetGenJetPDGID(int(genjet_hadf->at(ijet)));
    
    genjets.push_back(jet);
  }
  return genjets;
}


std::vector<KJet> SKTreeFiller::GetAllJets(){

  std::vector<KJet> jets;
  if(!SNUinput){

    for(std::vector<KJet>::iterator kit  = k_inputjets->begin(); kit != k_inputjets->end(); kit++){
      jets.push_back(*kit);
    }
    return jets;
  }

  for (UInt_t ijet=0; ijet< jet_eta->size(); ijet++) {
    KJet jet;
    if(jet_pt->at(ijet) != jet_pt->at(ijet)) continue;

    if(IsData){
      jet.SetPtEtaPhiE(jet_pt->at(ijet), jet_eta->at(ijet), jet_phi->at(ijet), jet_energy->at(ijet));
    }
    else{
      jet.SetPtEtaPhiE(jet_pt->at(ijet), jet_eta->at(ijet), jet_phi->at(ijet), jet_energy->at(ijet));

      // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution      
      /// Measurements show that the jet energy resolution (JER) in data is worse than in the simulation and the jets in MC need to be smeared to describe the data.

      jet*= jet_smearedRes->at(ijet); // JSKIM
      jet.SetIsMCSmeared(true);
    }
    jet.SetJetPassLooseID(jet_isLoose->at(ijet));
    jet.SetJetPassTightID(jet_isTight->at(ijet));
    jet.SetJetPassTightLepVetoID(jet_isTightLepVetoJetID->at(ijet));
    
    jet.SetJetPileupIDMVA(jet_PileupJetId->at(ijet));

    if(jet_PileupJetId){ 
      if(std::abs(jet_eta->at(ijet)) < 2.6){
	if(jet_PileupJetId->at(ijet) > 0.3) jet.SetJetPileupIDLooseWP(true);
	else jet.SetJetPileupIDLooseWP(false);
	if(jet_PileupJetId->at(ijet) > 0.7) jet.SetJetPileupIDMediumWP(true);
	else jet.SetJetPileupIDMediumWP(false);
	if(jet_PileupJetId->at(ijet) > 0.9)jet.SetJetPileupIDTightWP(true);
	else jet.SetJetPileupIDTightWP(false);
      }
      else{
	if(jet_PileupJetId->at(ijet) > -0.55) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(jet_PileupJetId->at(ijet) > -0.3) jet.SetJetPileupIDMediumWP(true);
	else jet.SetJetPileupIDMediumWP(false);
        if(jet_PileupJetId->at(ijet) > -0.1)jet.SetJetPileupIDTightWP(true);
	else jet.SetJetPileupIDTightWP(false);
      }
    }
    
    
    /// BTAG variables
    if(jet_CSVInclV2) jet.SetBTagInfo(snu::KJet::CSVv2, jet_CSVInclV2->at(ijet));
    if(jet_CMVAV2)    jet.SetBTagInfo(snu::KJet::cMVAv2, jet_CMVAV2->at(ijet));
    if(jet_JetProbBJet)  jet.SetBTagInfo(snu::KJet::JETPROB, jet_JetProbBJet->at(ijet)); 
    
    //if(jet_iCSVCvsL) {
    //      if(jet_iCSVCvsL->size() > 0)jet.SetCTagInfo(snu::KJet::iCSVCvsL, jet_iCSVCvsL->at(ijet));
    //    }
    if(jet_CCvsLT){
      if(jet_CCvsLT->size() > 0) jet.SetCTagInfo(snu::KJet::CCvsLT, jet_CCvsLT->at(ijet));
    }
    if(jet_CCvsBT){
      if(jet_CCvsBT->size() > 0)jet.SetCTagInfo(snu::KJet::CCvsBT, jet_CCvsBT->at(ijet));
    }
    jet.SetVtxMass(jet_vtxMass->at(ijet));
    jet.SetVtx3DVal(jet_vtx3DVal->at(ijet));
    jet.SetVtx3DSig(jet_vtx3DSig->at(ijet));
    jet.SetVtxNTracks(jet_vtxNtracks->at(ijet));
    
    // flavour
    jet.SetJetPartonFlavour(jet_partonFlavour->at(ijet));
    jet.SetJetHadronFlavour(jet_hadronFlavour->at(ijet));    
    jet.SetJetPartonPdgId(jet_partonPdgId->at(ijet));
    
    jet.SetJetChargedEmEF(jet_chargedEmEnergyFraction->at(ijet));
    /// JEC and uncertainties
    jet.SetJetScaledDownEnergy(jet_shiftedEnDown->at(ijet));
    jet.SetJetScaledUpEnergy(jet_shiftedEnUp->at(ijet));
    jet.SetSmearedResDown(jet_smearedResDown->at(ijet));
    jet.SetSmearedResUp(jet_smearedResUp->at(ijet));
    jet.SetSmearedRes(jet_smearedRes->at(ijet));
    
    if(jet_l1jetcorr){
      jet.SetJetRawPt(jet_rawpt->at(ijet));
      jet.SetJetRawEnergy(jet_rawenergy->at(ijet));
      jet.SetL1JetCorr(jet_l1jetcorr->at(ijet));
      jet.SetL2JetCorr(jet_l2jetcorr->at(ijet));
      jet.SetL3JetCorr(jet_l3jetcorr->at(ijet));
      jet.SetL2L3ResJetCorr(jet_l2l3resjetcorr->at(ijet));
      jet.SetJetArea(jet_area->at(ijet));
    }

    jets.push_back(jet);
  }// end of jet 
  
  
  std::sort( jets.begin(), jets.end(), isHigherPt );
  
  m_logger << DEBUG << "PFJet size = " << jets.size() << SNULogger::endmsg;
  return jets;
}



std::vector<KFatJet> SKTreeFiller::GetAllFatJets(){

  std::vector<KFatJet> fatjets;

  if(k_cat_version <  7) return fatjets;

  if(!SNUinput){

    for(std::vector<KFatJet>::iterator kit  = k_inputfatjets->begin(); kit != k_inputfatjets->end(); kit++){
      fatjets.push_back(*kit);
    }
    return fatjets;
  }

  for (UInt_t ijet=0; ijet< fatjet_eta->size(); ijet++) {
    KFatJet jet;
    if(fatjet_pt->at(ijet) != fatjet_pt->at(ijet)) continue;
    jet.SetPtEtaPhiE(fatjet_pt->at(ijet), fatjet_eta->at(ijet), fatjet_phi->at(ijet), fatjet_energy->at(ijet));

    jet.SetJetPassLooseID(fatjet_isLoose->at(ijet));
    jet.SetJetPassTightID(fatjet_isTight->at(ijet));
    jet.SetJetPassTightLepVetoID(fatjet_isTightLepVetoJetID->at(ijet));

    jet.SetJetPileupIDMVA(fatjet_PileupJetId->at(ijet));

    if(fatjet_PileupJetId){
      if(std::abs(fatjet_eta->at(ijet)) < 2.6){
        if(fatjet_PileupJetId->at(ijet) > 0.3) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(fatjet_PileupJetId->at(ijet) > 0.7) jet.SetJetPileupIDMediumWP(true);
        else jet.SetJetPileupIDMediumWP(false);
        if(fatjet_PileupJetId->at(ijet) > 0.9)jet.SetJetPileupIDTightWP(true);
        else jet.SetJetPileupIDTightWP(false);
      }
      else{
        if(fatjet_PileupJetId->at(ijet) > -0.55) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(fatjet_PileupJetId->at(ijet) > -0.3) jet.SetJetPileupIDMediumWP(true);
        else jet.SetJetPileupIDMediumWP(false);
        if(fatjet_PileupJetId->at(ijet) > -0.1)jet.SetJetPileupIDTightWP(true);
        else jet.SetJetPileupIDTightWP(false);
      }
    }


    /// BTAG variables                                                                                                                                                                                                                                                                                          
    if(fatjet_CSVInclV2) jet.SetBTagInfo(snu::KFatJet::CSVv2, fatjet_CSVInclV2->at(ijet));
    if(fatjet_CMVAV2)    jet.SetBTagInfo(snu::KFatJet::cMVAv2, fatjet_CMVAV2->at(ijet));
    if(fatjet_JetProbBJet)  jet.SetBTagInfo(snu::KFatJet::JETPROB, fatjet_JetProbBJet->at(ijet));

    if(fatjet_CCvsLT){
      if(fatjet_CCvsLT->size() > 0) jet.SetCTagInfo(snu::KFatJet::CCvsLT, fatjet_CCvsLT->at(ijet));
    }
    if(fatjet_CCvsBT){
      if(fatjet_CCvsBT->size() > 0)jet.SetCTagInfo(snu::KFatJet::CCvsBT, fatjet_CCvsBT->at(ijet));
    }
    jet.SetVtxMass(fatjet_vtxMass->at(ijet));
    jet.SetVtx3DVal(fatjet_vtx3DVal->at(ijet));
    jet.SetVtx3DSig(fatjet_vtx3DSig->at(ijet));
    jet.SetVtxNTracks(fatjet_vtxNtracks->at(ijet));

    // flavour                                                                                                                                                                                                                                                                                                  
    jet.SetJetPartonFlavour(fatjet_partonFlavour->at(ijet));
    jet.SetJetHadronFlavour(fatjet_hadronFlavour->at(ijet));
    jet.SetJetPartonPdgId(fatjet_partonPdgId->at(ijet));

    jet.SetJetChargedEmEF(fatjet_chargedEmEnergyFraction->at(ijet));

    jet.SetJetScaledDownEnergy(fatjet_shiftedEnDown->at(ijet));
    jet.SetJetScaledUpEnergy(fatjet_shiftedEnUp->at(ijet));
    jet.SetSmearedResDown(fatjet_smearedResDown->at(ijet));
    jet.SetSmearedResUp(fatjet_smearedResUp->at(ijet));
    jet.SetSmearedRes(fatjet_smearedRes->at(ijet));


    jet.SetTau1(fatjet_tau1->at(ijet));
    jet.SetTau2(fatjet_tau2->at(ijet));
    jet.SetTau3(fatjet_tau3->at(ijet));

    jet.SetPrunedMass(fatjet_prunedmass->at(ijet));
    jet.SetSoftDropMass(fatjet_softdropmass->at(ijet));

    jet.SetPuppiTau1(fatjet_puppi_tau1->at(ijet));
    jet.SetPuppiTau2(fatjet_puppi_tau2->at(ijet));
    jet.SetPuppiTau3(fatjet_puppi_tau3->at(ijet));
    jet.SetPuppiPt(fatjet_puppi_pt->at(ijet));
    jet.SetPuppiEta(fatjet_puppi_eta->at(ijet));
    jet.SetPuppiPhi(fatjet_puppi_phi->at(ijet));
    jet.SetPuppiM(fatjet_puppi_m->at(ijet));
    if(fatjet_l1jetcorr){
      jet.SetL1JetCorr(fatjet_l1jetcorr->at(ijet));
      jet.SetL2JetCorr(fatjet_l2jetcorr->at(ijet));
      jet.SetL3JetCorr(fatjet_l3jetcorr->at(ijet));
      jet.SetL2L3ResJetCorr(fatjet_l2l3resjetcorr->at(ijet));
      jet.SetJetArea(fatjet_area->at(ijet));
    }
      
    fatjets.push_back(jet);
  }// end of jet                                                   
  std::sort( fatjets.begin(), fatjets.end(), isHigherPt );

  m_logger << DEBUG << "PFJet size = " << fatjets.size() << SNULogger::endmsg;
  return fatjets;
}



std::vector<KMuon> SKTreeFiller::GetAllMuons(){

  std::vector<KMuon> muons ;
  
  if(!SNUinput){
    for(std::vector<KMuon>::iterator kit  = k_inputmuons->begin(); kit != k_inputmuons->end(); kit++){
      muons.push_back(*kit);
    }  
    return muons;
  }

  m_logger << DEBUG << "Filling Muons" << SNULogger::endmsg;

  vector<int> matched_truth;
  for (UInt_t ilep=0; ilep< muon_eta->size(); ilep++) {
    KMuon muon;
    if(muon_pt->at(ilep) != muon_pt->at(ilep)) continue;
    m_logger << DEBUG << "Filling global pt/eta ... " << SNULogger::endmsg;
   
    muon.SetTrigMatch(muon_trigmatch->at(ilep));
      
    /// GENERAL
    
    muon.SetISPF(muon_isPF->at(ilep));
    muon.SetIsGlobal(muon_isGlobal->at(ilep));

    muon.SetIsTracker(muon_isTracker->at(ilep));
    muon.SetIsLoose(muon_isLoose->at(ilep));
    muon.SetIsMedium(muon_isMedium->at(ilep));
    muon.SetIsTight(muon_isTight->at(ilep));
    muon.SetIsSoft(muon_isSoft->at(ilep));

    if(muon_shiftedEup){
      muon.SetShiftedEUp(muon_shiftedEup->at(ilep));
      muon.SetShiftedEDown(muon_shiftedEdown->at(ilep));
    }
    
    
    muon.SetPtEtaPhiE(muon_pt->at(ilep), muon_eta->at(ilep),muon_phi->at(ilep), muon_energy->at(ilep));
    if(k_cat_version > 4){
      muon.SetRochEta(muon_roch_eta->at(ilep));
      muon.SetRochPhi(muon_roch_phi->at(ilep));
      muon.SetRochE(muon_roch_energy->at(ilep));
      muon.SetRochM(muon_roch_m->at(ilep));
    }
    else{
      muon.SetRochPt(muon_pt->at(ilep));
      muon.SetRochEta(muon_eta->at(ilep));
      muon.SetRochPhi(muon_phi->at(ilep));
      muon.SetRochE(muon_energy->at(ilep));
      muon.SetRochM(muon_m->at(ilep));
    }
    muon.SetCharge(muon_q->at(ilep));
     
    m_logger << DEBUG << "Filling ms pt/eta ... " << SNULogger::endmsg;
 
    muon.SetRelIso(0.3,muon_relIso03->at(ilep));
    muon.SetRelIso(0.4,muon_relIso04->at(ilep));
    if(muon_minirelIso)muon.SetRelMiniIso(muon_minirelIso->at(ilep));

    if(k_cat_version  > 7){
      muon.SetMiniAODPt(muon_pt->at(ilep));
      muon.SetMiniAODRelIso(0.3,muon_relIso03->at(ilep));
      muon.SetMiniAODRelIso(0.4,muon_relIso04->at(ilep));
      muon.SetIsRochesterCorrected(false);
    }
    muon.Setdz(muon_dz->at(ilep));
    muon.Setdxy(muon_dxy->at(ilep));
    if(muon_sigdxy)muon.Setdxy_sig(muon_sigdxy->at(ilep));
    //// chi2
    muon.SetGlobalchi2( muon_normchi->at(ilep));
        
    /// hits
    muon.SetValidHits( muon_validhits->at(ilep));
    muon.SetPixelValidHits( muon_validpixhits->at(ilep));
    muon.SetValidStations( muon_matchedstations->at(ilep));
    muon.SetLayersWithMeasurement ( muon_trackerlayers->at(ilep));
    
    muon.SetMCMatched(muon_matched->at(ilep));


    muon.SetTrackVx(muon_x->at(ilep));
    muon.SetTrackVy(muon_y->at(ilep));
    muon.SetTrackVz(muon_z->at(ilep));

    //// Set Is ChargeFlip
    bool isprompt= false;
    bool from_tau = false;

    int mother_index=-1;
    int mother_pdgid=-1;
    int matched_index=-1;
    int mc_pdgid=-1;
    bool matched_in_Dr=false;
    
    int          mutype=0;

    if(k_cat_version > 3){

    if(gen_pt){
      float min_Dr=0.1;

      for (UInt_t it=0; it< gen_pt->size(); it++ ){
        if(gen_mother_index->at(it) <= 0)continue;
        if(gen_mother_index->at(it) >= int(gen_pt->size()))continue;
	
	if(gen_pt->at(it) < 5.) continue;

	double match_pt =muon_pt->at(ilep);
	double match_eta =muon_eta->at(ilep);
	double match_phi =muon_phi->at(ilep);
	double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));
	/// Matching using instructions on
        /// https://indico.cern.ch/event/292928/contributions/1650088/attachments/547844/755123/talk_electron_contribution.pdf
        ///

	bool already_matched=false;
	for(unsigned int im=0; im < matched_truth.size();im++){
	  if(it == unsigned(matched_truth.at(im))) already_matched=true;
	}

	///	if(already_matched) continue;

        /// Match requires to status 1 muon
        if(gen_status->at(it) != 1) continue;
        if(fabs(gen_PID->at(it)) != 13) continue;
	
	if(matched_in_Dr){
	  if(gen_mother_index->at(it) == gen_mother_index->at(matched_index)){
	    if (dr < 0.1){
	      // in case 2+ electrons from same mother electron (conversion) also match in pt
	      if( fabs(gen_pt->at(it)-muon_pt->at(ilep)) < fabs(gen_pt->at(matched_index)-muon_pt->at(ilep))) matched_index=it;
	    }
	  }
	  else   if (dr < min_Dr){
	    /// find closest match in dR to status 1
	    if( (fabs(match_pt - gen_pt->at(it))/gen_pt->at(it)) < 2.) {
	      matched_in_Dr=true;
	      min_Dr= dr;
	      /// set index of matched status 1 muon
	      matched_index=it;
	    }
	  }
	}
	else if (dr < min_Dr){
	  if( (fabs(match_pt - gen_pt->at(it))/gen_pt->at(it)) < 2.) {
	    /// find closest match in dR to status 1
	    matched_in_Dr=true;
	    min_Dr= dr;
	    
	    /// set index of matched status 1 muon
	    matched_index=it;
	  }
	}
      }// end of gen loop to find status 1 muon

      ///// treat case where there is a matched status 1 muon:
      //// classify into prompt:Fake:FromTau

      
      if(matched_in_Dr){
        /// Find closest non muon ancesteror
        float pdgid = gen_PID->at(matched_index);

        // mc_pdgid = closest matched status 1 pdgid
        mc_pdgid= int(pdgid);

        // mindex = mother index: will loop to find first non el mother
        int mindex= matched_index;

        while ( (fabs(gen_PID->at(mindex)) == 13)) {
          pdgid = gen_PID->at(mindex);
          mindex=gen_mother_index->at(mindex);
        }
        /// pdgid is now of muon from non muon mother
        //  mindex = index for mother of non muon ancestor

        if( (fabs(gen_PID->at(mindex)) == 23) || (fabs(gen_PID->at(mindex)) == 24)) {
	  /// Check if it is a chargeflip.
          mutype=1;
	  
          if(pdgid * muon_q->at(ilep) > 0 )     muon.SetIsChargeFlip(true);
          else     muon.SetIsChargeFlip(false);

	  int n_mu_from_mother=0;
	  for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	    if(gen_mother_index->at(itx) <= 0)continue;
	    if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(itx) < 0.001) continue;
	    if(fabs(gen_PID->at(itx)) ==13) {
	      if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) n_mu_from_mother++;
	    }
	  }
	  if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
	  if(n_mu_from_mother == 5)  muon.SetIsPhotonConversion(true);
          mother_index=mindex;
          mother_pdgid=gen_PID->at(mindex);
          isprompt=true; /// means is prompt
        }/// end of Z/W
        else {
          if(gen_status->at(mindex) == 2){
            if(fabs(gen_PID->at(mindex)) > 50) {isprompt=false; mother_pdgid=gen_PID->at(mindex); mother_index=mindex; from_tau=false;
	      mutype=2;
	    }
	    else {
	      isprompt=true;
	      mutype=3;

	      if(fabs(gen_PID->at(mindex)) == 22){
		if(fabs(gen_PID->at(gen_mother_index->at(mindex))) > 50){
                  mutype=4;
                }
		else                   mutype=5;
              } 
	    }
	    
            if(fabs(gen_PID->at(mindex)) == 15){
              isprompt=true; mother_pdgid=gen_PID->at(mindex);  mother_index=mindex; from_tau=true;
              // Check if el from tau  is CF
	      mutype=6;
	      if(pdgid * muon_q->at(ilep) > 0 )     muon.SetIsChargeFlip(true);
              else     muon.SetIsChargeFlip(false);

	      int n_mu_from_mother=0;
	      for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		if(gen_mother_index->at(itx) <= 0)continue;
		if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
		if(gen_pt->at(itx) < 0.001) continue;
		if(fabs(gen_PID->at(itx)) ==13) {
		  if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) n_mu_from_mother++;
		}
	      }
	      if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
	      if(n_mu_from_mother == 5)  muon.SetIsPhotonConversion(true);
	     
	      if(fabs(gen_PID->at(gen_mother_index->at(mother_index))) > 50) {isprompt=false; mutype=7;}
            }
          }/// end of status 2 check
          else {
            /// using new method for matching: These events are set as prompt
            isprompt=true;mother_pdgid=-99999; mother_index=mindex; from_tau=false;
	    mutype=8;
            if(pdgid * muon_q->at(ilep) > 0 )    muon.SetIsChargeFlip(true);
            else     muon.SetIsChargeFlip(false);
	    int n_mu_from_mother=0;
            vector<KTruth> vmu_tmp;
            bool neutrino_invertex(false);

	    for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	      if(gen_mother_index->at(itx) <= 0)continue;
	      if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
	      if(gen_pt->at(itx) < 0.001) continue;
	      if(fabs(gen_PID->at(itx)) ==13) {
		if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) n_mu_from_mother++;
		if(gen_status->at(itx) ==1){
		  KTruth truthmu;
		  truthmu.SetPtEtaPhiE(gen_pt->at(itx), gen_eta->at(itx), gen_phi->at(itx), gen_energy->at(itx));
		  vmu_tmp.push_back(truthmu);
		}
	      }
	      if(fabs(gen_PID->at(itx)) ==14) {
		
		int index_mother_nu=gen_mother_index->at(itx);
		while (fabs(index_mother_nu) ==14){
		  index_mother_nu=gen_mother_index->at(index_mother_nu);
		}
		
		if(index_mother_nu == mindex) {
		  neutrino_invertex=true;
		}
	      }
	    } // end of truth loop to check Conv                                                                                                                                                                                                  
	    
	    if(neutrino_invertex) mutype=9;
	    
	    if(vmu_tmp.size() ==2) {
	      KParticle ll = vmu_tmp[0] + vmu_tmp[1];
	      if(fabs(ll.M()) < 5.) mutype=10;
	    }
	    
	    
	    
	    if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
	    if(n_mu_from_mother == 5)  muon.SetIsPhotonConversion(true);

	    /// speacial treatment for signal                                                                                                                                                                                                                                   

            if( fabs(gen_PID->at(mindex))>= 9900012 &&  fabs(gen_PID->at(mindex) < 9900025 ))mother_pdgid= gen_PID->at(mindex);


          }
        }
      }      /// In case no status 1 muon is found : classify muon fake
      else{
	if(muon_q->size() ==2){
	  for (UInt_t itxx=0; itxx< gen_pt->size(); itxx++ ){
            if(gen_mother_index->at(itxx) <= 0)continue;
            if(gen_mother_index->at(itxx) >= int(gen_pt->size()))continue;
            if(gen_pt->at(itxx) < 0.001) continue;
	    //cout << itxx << " " << gen_PID->at(itxx) << " " << gen_PID->at(gen_mother_index->at(itxx)) << " " << gen_eta->at(itxx) << " " << gen_phi->at(itxx) << endl;
	  }
	}
        if(gen_pt){
	  for (UInt_t it=0; it< gen_pt->size(); it++ ){
            if(gen_mother_index->at(it) <= 0)continue;
            if(gen_mother_index->at(it) >= int(gen_pt->size()))continue;
            if(gen_pt->at(it) < 0.001) continue;

	    bool already_matched=false;
	    for(unsigned int im=0; im < matched_truth.size();im++){
	      if(it == unsigned(matched_truth.at(im))) already_matched=true;
	    }
	    //if(already_matched) continue;

            double match_eta =muon_eta->at(ilep);
            double match_phi =muon_phi->at(ilep);
            double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));

            if (dr <0.1){
              matched_in_Dr=true;
              int mindex= gen_mother_index->at(it);
              float pdgid = gen_PID->at(it);
	      
	      /// Unlikely to have mother as muon but just in case
              while ( (fabs(gen_PID->at(mindex)) == 13)) {
                mindex=gen_mother_index->at(mindex);
              }
              // isprompt = false since it failed status 1 matching
              isprompt=false;
              /// mother index of first non muon
              mother_pdgid=gen_PID->at(mindex);
              mother_index=mindex;
              matched_index = it;
	      mc_pdgid= int(pdgid);
              if(fabs(pdgid) == 22) {
		if(fabs(mother_pdgid) > 50) mutype=11;
                else mutype=12;
		
		from_tau=false;
                break;
              }

              int grandmother =0;
              if(gen_mother_index->at(mindex) > 0)grandmother=gen_PID->at(gen_mother_index->at(mindex));

              if(fabs(pdgid) == 15){
		from_tau=true;
		if(fabs(mother_pdgid) > 50) mutype=13;
		else mutype=14;
		if(fabs(mother_pdgid) ==15){
		  if(fabs(grandmother) > 50) mutype=-13;
		  else mutype=-14;

		}
	      }
	      if(fabs(gen_PID->at(gen_mother_index->at(it))) == 13){
		mutype=15;
		int n_mu_from_mother=0;
		for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		  if(gen_mother_index->at(itx) <= 0)continue;
		  if(gen_mother_index->at(itx) >= int(gen_pt->size()))continue;
		  if(gen_pt->at(itx) < 0.001) continue;
		  if(fabs(gen_PID->at(itx)) ==13) {
		    if(gen_mother_index->at(itx) == gen_mother_index->at(matched_index)) n_mu_from_mother++;
		  }
		}
		if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
		if(n_mu_from_mother == 5)  muon.SetIsPhotonConversion(true);
	      }else{
		muon.SetIsPhotonConversion(false);
	      }
	      
	      if(fabs(pdgid) == 1){
                mutype=16;
              }
              else if(fabs(pdgid) == 2){
                mutype=17;
              }
              else if(fabs(pdgid) == 3){
               mutype=18;
              }
              else if(fabs(pdgid) == 4){
                mutype=19;
              }
              else if(fabs(pdgid) == 5){
                mutype=20;
              }
              else if(fabs(pdgid) == 21){
                mutype=21;
              }
	      
	      else if(fabs(pdgid) == 211){
		mutype=22;
	      }
	      else if(fabs(pdgid) == 310){
		mutype=23;
	      }
	      else if(fabs(pdgid) == 431){
		mutype=24;
	    }
	      else if(fabs(pdgid) == 13){
		mutype=25;
	      }
	      else if(fabs(pdgid) == 311){
                mutype=26;
              }
	      else if(fabs(pdgid) > 500 && fabs(pdgid) < 600){
                mutype=27;
              }
	      else mutype=28;

	    }// dr req
          }// loop over gen vector
        }// require gen info
      }// no status 1 match
    }
    
    /// matched_index is index which matches reco muon with smallest dR
    /// - If multiple status 1 muons are matched look at closest in pt
    /// - In no status 1 is matched set as not prompt but look for closest particle in dR 
    /// - In no particles within dR < 0.1 matched_in_Dr = false
    matched_truth.push_back(matched_index);


    muon.SetType(mutype);

    if(!matched_in_Dr){
      muon.SetMCMatched(false);
      muon.SetIsFromTau(false);
      muon.SetMCMatchedPdgId(-1);
      muon.SetMotherPdgId(0);
      muon.SetMotherTruthIndex(-1);
      muon.SetMCTruthIndex(-1);
    }
    else{

      muon.SetMCMatched(isprompt);
      muon.SetIsFromTau(from_tau);
      muon.SetMotherPdgId(mother_pdgid);
      muon.SetMCMatchedPdgId(mc_pdgid);
      muon.SetMotherTruthIndex(mother_index);
      muon.SetMCTruthIndex(matched_index);
      if(gen_status->at(matched_index)==1)muon.SetIsPromptFlag(gen_isprompt->at(matched_index));

    }
    }

    /// Fill vector
    muons.push_back(muon);
  }
  
  std::sort( muons.begin(), muons.end(), isHigherPt );
  m_logger << DEBUG << "End of Muon Filling" << SNULogger::endmsg;
  return muons;
}

  

std::vector<snu::KTruth>   SKTreeFiller::GetTruthParticles(int np){
  
  m_logger << DEBUG << "Filling Truth" << SNULogger::endmsg;
  std::vector<snu::KTruth> vtruth;

  if(IsData) return vtruth;

  int counter=0;

  if(!SNUinput){

    for(std::vector<KTruth>::iterator kit  = k_inputtruth->begin(); kit != k_inputtruth->end(); kit++, counter++){
      if(counter == np)  break;
      vtruth.push_back(*kit);
    }
    return vtruth;
  }
  
  m_logger << DEBUG << "Filling truth Info: " << gen_pt->size() << SNULogger::endmsg;

  for (UInt_t it=0; it< gen_pt->size(); it++ , counter++) {
    
    if(counter == np)  break;
    KTruth truthp;
    truthp.SetPtEtaPhiE(double(gen_pt->at(it)), double(gen_eta->at(it)), double(gen_phi->at(it)), double(gen_energy->at(it)));
    truthp.SetParticlePdgId(gen_mother_PID->at(it));
    truthp.SetParticleStatus(gen_status->at(it));
    truthp.SetParticleIndexMother(gen_mother_index->at(it));
    
    if(k_cat_version > 3){
      // To save space set a single int as the flag. 
      // 
      int truth_flag = 0;
      if(gen_isprompt->at(it)) truth_flag+=1;
      if(gen_isdecayedleptonhadron->at(it)) truth_flag+=10;
      if(gen_istaudecayproduct->at(it)) truth_flag+=100;
      if(gen_isprompttaudecayproduct->at(it)) truth_flag+=1000;
      if(gen_isdirecthadrondecayproduct->at(it)) truth_flag+=10000;
      if(gen_ishardprocess->at(it)) truth_flag+=100000;
      if(gen_fromhardprocess->at(it)) truth_flag+=1000000;
      if(gen_fromhardprocess_beforeFSR->at(it)) truth_flag+=10000000;
      truthp.SetStatusFlag(truth_flag);
    }
    
    vtruth.push_back(truthp);  
  }
  
  return vtruth;
}
