#ifndef _SKTree_KElectron_H__
#define _SKTree_KElectron_H__

/// Local includes
#include "KParticle.h"

#include <iosfwd>
#include <string>
#include "TLorentzVector.h"

namespace snu {
  
  class KElectron : public KParticle {
  public:


    enum ElectronType{PROMPT=0,
		      FAKE=1,
		      PHOTONFAKE=2,
		      CONV_CF=3,
		      CONV_NONECF=4,
		      CF=5,
    };
    enum ElectronMotherType{none=0,
			    Z=1,
			    W=2,
			    ZorW=3,
			    pion=4,
    };

    KElectron();
    
    ///Copy constructor
    KElectron(const KElectron& el);
    
    ///Destructor    
    virtual ~KElectron() ;

    KElectron& operator= (const KElectron& obj);
    
    
    float ScaleFactor(const std::string& name, int sign) const ;
    Bool_t PassTrigMVAHNLoose() const ;
    Bool_t PassTrigMVAHNTightvX(float mvacut) const;
    Bool_t PassTrigMVAHNTight() const ;
    Bool_t PassTrigMVAHNTightv2(float pt1, float pt2, float mvacut, float mvacut_ec) const ;
    Bool_t PassTrigMVAHNTightv3() const ;
    Bool_t PassTrigMVAHNTightv4() const ;
    Bool_t PassTrigMVAGENTTight() const ;
    Bool_t PassTrigMVAGENTLoose() const;

    // set kinematic variables
    void SetSCEta(Double_t sceta);

    void SetSmearFactor(Double_t smear);

    
    //##### NOTE charge/pt/eta/phi use tlv class
    
    /// Trigger matching
    
    /// MVA
    void SetMVA(double mva);
    void SetZZMVA(double zzmva);
    void SetGSF(double pt, double eta, double phi, double m);
    void SetSC(double eta, double phi);
    
    void SetElVariables(double electron_etaWidth, double electron_phiWidth, double electron_dEtaIn, double electron_dEtaInSeed, double electron_dPhiIn, double electron_sigmaIEtaIEta, double electron_Full5x5_SigmaIEtaIEta, double electron_HoverE, double electron_fbrem, double electron_eOverP);
    //// set   vertex variables

    void SetType(int eltype);

    void Setdxy(Double_t d_xy);
    void Setdxy_sig(Double_t d_xysig);
    void Setdz(Double_t d_z);
    
    void SetSNUID(int id);
    void SetPassVeto(Bool_t pass);
    void SetPassLoose(Bool_t pass);
    void SetPassMedium(Bool_t pass);
    void SetPassTight(Bool_t pass);
    void SetPassHLT(Bool_t pass);
    void SetPassHEEP(Bool_t pass);
    
    

    void SetPassMVATrigMedium(Bool_t pass);
    void SetPassMVATrigTight(Bool_t pass);
    void SetPassMVANoTrigMedium(Bool_t pass);
    void SetPassMVANoTrigTight(Bool_t pass);
    void SetPassMVAZZ(Bool_t pass);
    
    void SetIsPF(Bool_t ispf);
    void SetIsChargeFlip(Bool_t iscf);
    void SetIsPhotonConversion(Bool_t isconv);
    void SetIsFromTau(Bool_t istau);
    void SetIsMCMatched(Bool_t ismatch);
    
    void SetMCMatchedPdgId(Int_t pg);
    void SetMotherPdgId(Int_t pg);
    void SetMotherTruthIndex(Int_t mindex);
    void SetMCTruthIndex(Int_t t_index);
    /// set ISO variables
    void SetPFChargedHadronIso(Double_t cone,Double_t pf_ch);
    void SetPFPhotonIso(Double_t cone,Double_t pf_ph);
    void SetPFNeutralHadronIso(Double_t cone,Double_t pf_ne);

    
    void SetPFRelIsoRho(Double_t cone, Double_t pf_rel);
    void SetPFRelIsoBeta(Double_t cone, Double_t pf_rel);
    void SetPFRelMiniIso( Double_t pf_rel);
    void SetPFAbsIso(Double_t cone, Double_t pf_abs);

    
    // set charge variables
    void SetGsfCtfScPixCharge(bool gsfctfscpix_ch);

    /// set conversion variables
    void SetIsMCExternalConversion(Bool_t isconv);
    void SetHasMatchedConvPhot(Bool_t hasmatchConvPhot);
    void SetMissingHits(Int_t mhits);
    
    void SetScaleEUp(Double_t Eup);
    void SetScaleEDown(Double_t Edown);
    
    void SetSmearEUp(Double_t Eup);
    void SetSmearEDown(Double_t Edown);


    void SetTrkVx(Double_t trkvx);
    void SetTrkVy(Double_t trkvy);
    void SetTrkVz(Double_t trkvz);

    void SetIsoMVA(Double_t mva);
    void SetNonIsoMVA(Double_t mva);

    void SetIsPromptFlag(bool pflag);

    bool TriggerMatched(TString path);

    ///// Functions to call class variables
    
    inline Int_t GetType()  const {
      return k_eltype;
    }
    inline Bool_t IsPromptFlag() const {return k_isprompt;}

    inline Double_t SmearFactor() const {return k_smearfactor;}

    inline Double_t IsoMVA() const {return k_mva_iso;}
    inline Double_t NonIsoMVA() const {return k_mva_noniso;}


    inline Bool_t  IsEBFiducial() {return bool (fabs(SCEta()) < 1.442);}
    inline Bool_t  IsEB1() {return bool (fabs(SCEta()) < 0.8);}
    inline Bool_t  IsEB2() {return bool (fabs(SCEta()) < 1.479);}
    inline Bool_t  IsEE() {return bool (fabs(SCEta()) > 1.479 && fabs(SCEta()) < 2.50);}
    inline Bool_t  IsEEFiducial() {return bool (fabs(SCEta()) > 1.560 && fabs(SCEta()) < 2.50);}
      
    /// // Kinematic variables
    inline Double_t  SCEta() const {return k_sceta;}
    
    
    inline Int_t MissingHits() const {return k_missing_hits;}

    inline Double_t PtShiftedUp() const{ return k_pt_shifted_up;}
    inline Double_t PtShiftedDown() const{ return k_pt_shifted_down;}

    /// Trigger matching
    
    // ID variables
    inline Bool_t PassVeto() const{return pass_veto;}
    inline Bool_t PassLoose() const{return pass_loose;}
    inline Bool_t PassMedium() const{return pass_medium;}
    inline Bool_t PassTight() const{return pass_tight;}
    inline Bool_t PassHLTID() const{return pass_hltid;}

    // HEEP ID 6.0
    inline Bool_t PassHEEP() const{return pass_heep;}

    
    inline Bool_t MCMatched() const{
      if(k_is_fromtau) return true;
      return k_mc_matched;
    }

    inline Bool_t IsPF() const{return k_isPF;}
    inline Bool_t MCIsPrompt() const{return k_mc_matched;}
    inline Bool_t MCIsCF() const{return k_is_cf;}
    inline Bool_t MCIsFromConversion() const{return k_is_conv;}
    inline Bool_t MCIsExternalConversion() const{return k_in_conv;}
    inline Bool_t MCFromTau() const{return k_is_fromtau;}
    inline Int_t MCMatchedPdgId() const{return k_mc_pdgid;}
    inline Int_t MotherPdgId() const{return k_mother_pdgid;}
    inline Int_t MotherTruthIndex() const{return k_mother_index;}
    inline Int_t MCTruthIndex() const{return k_mc_index;}
    
    inline Int_t SNUID() const{return snu_id;}

    inline KElectron::ElectronType GetParticleType() const{ 
      if(k_is_conv&&k_is_cf) return KElectron::CONV_CF;
      if(k_is_conv&&!k_is_cf)   return KElectron::CONV_NONECF; 
      if(k_is_cf)  return KElectron::CF;
      if(k_mc_matched) return KElectron::PROMPT;
      if(k_mc_pdgid==22) return KElectron::PHOTONFAKE;
      return KElectron::FAKE;

    }

    inline KElectron::ElectronMotherType GetMotherType() const{
      if(k_mother_pdgid == 23) return KElectron::Z;
      if(fabs(k_mother_pdgid) == 24) return KElectron::W;
      if(k_mother_pdgid == -99999 ) return KElectron::ZorW;
      return  KElectron::pion;

    }
    // charge variables
    
    inline Bool_t GsfCtfScPixChargeConsistency()  const {return k_gsf_ctscpix_charge;}
    inline Bool_t GsfCtfScPixChargeConsistency()  const {return k_gsf_ctscpix_charge;}
    inline Bool_t GsfCtfScPixChargeConsistency()  const {return k_gsf_ctscpix_charge;}
    
    // Conversion variables
    inline Bool_t PassesConvVeto() const {return k_hasmatchconvphot;}
    

    
    // Isolation Variables
    inline Double_t PFChargedHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_chargedhad_iso03;
      else  if(cone == 0.4)  return k_pf_chargedhad_iso04;
      else return -999.;
    }
    
    inline Double_t PFPhotonIso(double cone) const {
      if(cone == 0.3)   return k_pf_photon_iso03;
      else  if(cone == 0.4)  return k_pf_photon_iso04;
      else return -999.;
    }
    
    inline Double_t PFNeutralHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_neutral_iso03;
      else  if(cone == 0.4)  return k_pf_neutral_iso04;
      else return -999.;
    }
    
		
    inline Double_t PTCone(double cone, double iso=0.08){
      float ptcone= this->Pt() * (1.+ std::max(0., k_rel_iso03 - iso));
      if(cone == 0.4)   ptcone= this->Pt() * (1.+ std::max(0., k_rel_iso04 - iso));
      if(ptcone< 10.) return -1.;
      return ptcone;
    }


    inline Double_t PFRelIso(double cone) const {
      if(cone == 0.3)   return k_rel_iso03;
      else  if(cone == 0.4)   return k_rel_iso04;
      else return -999.;
    }

    inline Double_t PFRelMiniIso() const { return k_rel_miniiso; }


    inline Double_t PFAbsIso(double cone) const {
      if(cone == 0.3)   return k_abs_iso03;
      else  if(cone == 0.4)   return k_abs_iso04;
      else return -999.;
    }

    /// VtxDist with vertex chosen to be primary   
    inline Double_t  dxy() const {return  k_dxy;}
    inline Double_t  dxySig() const {return  k_dxy_sig;}
    inline Double_t  dz() const {return  k_dz;}
    

    inline Double_t  TrkVx() const {return  k_trkvx;}
    inline Double_t  TrkVy() const {return  k_trkvy;}
    inline Double_t  TrkVz() const {return  k_trkvz;}

    inline TString TrigMatch() const{return k_trig_match;}

  protected:
    /// Reset function.                                                                  
    virtual void Reset();    
    
  private:
    /// decalre private functions

    Double_t k_pf_chargedhad_iso03, k_pf_photon_iso03, k_pf_neutral_iso03, k_pf_chargedhad_iso04, k_pf_photon_iso04, k_pf_neutral_iso04, k_rel_iso03, k_rel_iso04, k_rel_miniiso;
    Double_t k_abs_iso03, k_abs_iso04;
    Double_t k_electron_dxyVTX,k_electron_dzVTX, k_electron_dxy, k_electron_sigdxy, k_electron_dz, k_electron_ip3D, k_electron_sigip3D, k_electron_dxyBS, k_electron_dzBS;
    Double_t k_electron_AEff03, k_electron_chIso03, k_electron_nhIso03, k_electron_phIso03, k_electron_pcIso03, k_electron_puChIso03, k_electron_chIso04, k_electron_nhIso04, k_electron_phIso04, k_electron_pcIso04, k_electron_puChIso04;
    Double_t k_electron_relIsoCom03, k_electron_relIsoCom04, k_electron_relIsoBeta03, k_electron_relIsoBeta04, k_electron_relIsoRho03;

    
    Bool_t k_electron_isGsfCtfScPixChargeConsistent, electron_isGsfScPixChargeConsistent,electron_isGsfCtfChargeConsistent;
    Bool_t pass_hltid,pass_tight, pass_veto, pass_medium, pass_loose, k_mc_matched,  k_is_cf,k_is_conv, k_is_fromtau,k_isPF,k_hasmatchconvphot, pass_heep;
    
    Double_t k_pt_scale_up, k_pt_scale_down, k_pt_smear_up, k_pt_smear_down;
    Int_t snu_id,k_mother_pdgid, k_mc_pdgid,k_mother_index, k_mc_index;
    TString k_trig_match;
    Int_t k_eltype;
    
    Double_t k_mva_iso, k_mva_noniso;
    Int_t k_electron_mHits;
    Double_t  k_smear_up, k_smear_down,k_scale_up, k_scale_down;
    Double_t k_gsf_pt, k_gsf_eta, k_gsf_phi, k_gsf_charge;
    Double_t k_sceta, k_scphi;
    Double_t k_electron_dEtaIn, k_electron_dEtaInSeed, k_electron_dPhiIn, k_electron_sigmaIEtaIEta, k_electron_Full5x5_SigmaIEtaIEta,electron_HoverE,electron_fbrem,electron_eOverP, electron_energyEC, electron_Pnorm,electron_InvEminusInvP;
    Bool_t k_is_conv;
    Bool_t k_isprompt;


    ClassDef(KElectron,32);
  }; 
  
}//namespace snu

#endif
