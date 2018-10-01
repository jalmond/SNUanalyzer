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
    
    

    // set kinematic variables
    void SetSCEta(Double_t sceta);
    void SetSCPhi(Double_t scphi);


    
    //##### NOTE charge/pt/eta/phi use tlv class
    
    /// MVA
    void SetMVAIso(double mva);
    void SetMVANonIso(double zzmva);

    //    void SetGSF(double pt, double eta, double phi, double m);
    


    void SetType(int eltype);

    void SetIP2D(Double_t d_ip2d);
    void SetIP3D(Double_t d_ip3d);
    void SetSIP3D(Double_t d_sip3d);

    void Setdz(Double_t d_z);
    
    void SetPassVeto(Bool_t pass);
    void SetPassLoose(Bool_t pass);
    void SetPassMedium(Bool_t pass);
    void SetPassTight(Bool_t pass);
    void SetPassHEEP(Bool_t pass);
    
    
    void SetIsChargeFlip(Bool_t iscf);
    void SetIsPhotonConversion(Bool_t isconv);

    void SetIsFromTau(Bool_t istau);
    void SetIsMCMatched(Bool_t ismatch);
       void SetMCMatchedPdgId(Int_t pg);	void SetMotherPdgId(Int_t pg);
    void SetMotherTruthIndex(Int_t mindex);
    void SetMCTruthIndex(Int_t t_index);

    /// set ISO variables
    void SetPFChargedHadronIso(Double_t cone,Double_t pf_ch);
    void SetPFPhotonIso(Double_t cone,Double_t pf_ph);
    void SetPFNeutralHadronIso(Double_t cone,Double_t pf_ne);
    void SetPFRelIsoRho(Double_t cone, Double_t pf_rel);
    void SetPFRelIsoBeta(Double_t cone, Double_t pf_rel);
    void SetPFRelMiniIso( Double_t pf_rel);
  
    
    // set charge variables
    void SetGsfCtfScPixCharge(bool gsfctfscpix_ch);
    void SetGsfScPixCharge(bool gsfscpix_ch);
    void SetGsfCtfCharge(bool gsfctf_ch);


    /// set conversion variables
    void SetIsMCExternalConversion(Bool_t isconv);
    void SetHasMatchedConvPhot(Bool_t hasmatchConvPhot);
    void SetMissingHits(Int_t mhits);
    void SetEcalDriven(Int_t edriven);
    
    void SetEnUncorr(Double_t E);

    void SetScaleEUp(Double_t Eup);
    void SetScaleEDown(Double_t Edown);
    
    void SetSmearEUp(Double_t Eup);
    void SetSmearEDown(Double_t Edown);


    void SetIsPromptFlag(bool pflag);

    void SetElIDVariables(double  electron_Full5x5_SigmaIEtaIEta , double electron_dEtaInSeed ,double  k_electron_dPhiIn,double  k_electron_HoverE  ,double  k_electron_InvEminusInvP);

    ///// Functions to call class variables
    
    inline Int_t GetType()  const {
      return k_eltype;
    }
    inline Bool_t IsPromptFlag() const {return k_isprompt;}

    inline Double_t IsoMVA() const {return k_mva_iso;}
    inline Double_t NonIsoMVA() const {return k_mva_noniso;}


    inline Bool_t  IsEBFiducial() {return bool (fabs(SCEta()) < 1.442);}
    inline Bool_t  IsEB1() {return bool (fabs(SCEta()) < 0.8);}
    inline Bool_t  IsEB2() {return bool (fabs(SCEta()) < 1.479);}
    inline Bool_t  IsEE() {return bool (fabs(SCEta()) > 1.479 && fabs(SCEta()) < 2.50);}
    inline Bool_t  IsEEFiducial() {return bool (fabs(SCEta()) > 1.560 && fabs(SCEta()) < 2.50);}
      
    /// // Kinematic variables
    inline Double_t  SCEta() const {return k_sceta;}
    inline Double_t  SCPhi() const {return k_scphi;}
    
    
    inline Int_t MissingHits() const {return k_missing_hits;}
    inline Int_t EcalDriven() const {return k_ecaldriven;}

    inline Double_t PtScaleUp() const{ return k_pt_scale_up;}
    inline Double_t PtScaleDown() const{ return k_pt_scale_down;}
    inline Double_t PtSmearUp() const{ return k_pt_smear_up;}
    inline Double_t PtSmearDown() const{ return k_pt_smear_down;}

    
    // ID variables
    inline Bool_t PassVeto() const{return pass_veto;}
    inline Bool_t PassLoose() const{return pass_loose;}
    inline Bool_t PassMedium() const{return pass_medium;}
    inline Bool_t PassTight() const{return pass_tight;}

    // HEEP ID
    inline Bool_t PassHEEP() const{return pass_heep;}

    
    inline Bool_t MCMatched() const{
      if(k_is_fromtau) return true;
      return k_mc_matched;
    }

    inline Bool_t MCIsPrompt() const{return k_mc_matched;}
    inline Bool_t MCIsCF() const{return k_is_cf;}
    inline Bool_t MCIsFromConversion() const{return k_is_conv;}
    inline Bool_t MCIsExternalConversion() const{return k_is_conv;}
    inline Bool_t MCFromTau() const{return k_is_fromtau;}
    inline Int_t MCMatchedPdgId() const{return k_mc_pdgid;}
    inline Int_t MotherPdgId() const{return k_mother_pdgid;}
    inline Int_t MotherTruthIndex() const{return k_mother_index;}
    inline Int_t MCTruthIndex() const{return k_mc_index;}

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
    inline Bool_t GsfScPixChargeConsistency()  const {return k_gsf_scpix_charge;}
    inline Bool_t GsfCtfChargeConsistency()  const {return k_gsf_ct_charge;}

    
    // Conversion variables
    inline Bool_t PassesConvVeto() const {return k_hasmatchconvphot;}
    

    
    // Isolation Variables
    inline Double_t PFChargedHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_chargedhad_iso03;
      else return -999.;
    }
    
    inline Double_t PFPhotonIso(double cone) const {
      if(cone == 0.3)   return k_pf_photon_iso03;
      else return -999.;
    }
    
    inline Double_t PFNeutralHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_neutral_iso03;
      else return -999.;
    }
    
		
    inline Double_t PTCone(double iso=0.08){
      float ptcone= this->Pt() * (1.+ std::max(0., k_rel_iso03 - iso));

      if(ptcone< 10.) return -1.;
      return ptcone;
    }


    inline Double_t PFRelMiniIso() const { return k_rel_miniiso; }


    inline Double_t PFRelIsoBeta(double cone ) const {
      if(cone == 0.3)   return k_electron_relIsoBeta03;
      else return -999.;
    }
    inline Double_t PFRelIsoRho(double cone ) const {
      if(cone == 0.3)   return k_electron_relIsoRho03;
      else return -999.;
    }

    inline Double_t PFRelIso(double cone ) const {
      return  PFRelIsoRho(cone);

    }

    inline Double_t PFAbsIsoBeta(double cone) const {
      if(cone == 0.3)   return k_electron_relIsoBeta03 * this->Pt();
      else return -999.;
    }

    
    inline Double_t PFAbsIsoRho(double cone) const {
      if(cone == 0.3)   return k_electron_relIsoRho03 * this->Pt();
      else return -999.;
    }

    /// VtxDist with vertex chosen to be primary   
    inline Double_t  IP2D() const {return  k_ip2D;}
    inline Double_t  IP3D() const {return  k_ip3D;}
    inline Double_t  SIP3D() const {return  k_sip3D;}

    inline Double_t  dz() const {return  k_dz;}
    

    inline Double_t  MVAIso() const {return  k_mva_iso;}
    inline Double_t  MVANonIso() const {return  k_mva_noniso;}


    inline Double_t  Full5x5_SigmaIEtaIEta() const {return  k_electron_Full5x5_SigmaIEtaIEta;}
    inline Double_t  dEtaInSeed() const {return  k_electron_dEtaInSeed;}
    inline Double_t  dPhiIn() const {return  k_electron_dPhiIn;}
    inline Double_t  HoverE() const {return  k_electron_HoverE;}
    inline Double_t  InvEminusInvP() const {return  k_electron_InvEminusInvP;}



  protected:
    /// Reset function.                                                                  
    virtual void Reset();    
    
  private:
    /// decalre private functions

    Double_t k_pf_chargedhad_iso03, k_pf_photon_iso03, k_pf_neutral_iso03, k_rel_iso03, k_rel_miniiso;

    Double_t k_ip2D, k_ip3D, k_sip3D, k_dz;


    Double_t  k_electron_relIsoBeta03, k_electron_relIsoRho03;

    
    Bool_t k_gsf_ctscpix_charge, k_gsf_scpix_charge, k_gsf_ct_charge;

    Bool_t pass_hltid,pass_tight, pass_veto, pass_medium, pass_loose, k_mc_matched,  k_is_cf,k_is_conv, k_is_fromtau,k_hasmatchconvphot, pass_heep;
    
    Double_t k_pt_scale_up, k_pt_scale_down, k_pt_smear_up, k_pt_smear_down;
  
  Int_t snu_id,k_mother_pdgid, k_mc_pdgid,k_mother_index, k_mc_index;

    Int_t k_eltype;
    
    Double_t k_mva_iso, k_mva_noniso;

    Int_t k_missing_hits, k_ecaldriven;

    Double_t  k_smear_up, k_smear_down,k_scale_up, k_scale_down,k_en_uncorr;
 
    //    Double_t k_gsf_pt, k_gsf_eta, k_gsf_phi, k_gsf_charge;

    Double_t k_sceta, k_scphi;


    Bool_t k_isprompt;


    Double_t k_electron_Full5x5_SigmaIEtaIEta, k_electron_dEtaInSeed, k_electron_dPhiIn , k_electron_HoverE, k_electron_InvEminusInvP;

    ClassDef(KElectron,33);
  }; 
  
}//namespace snu

#endif
