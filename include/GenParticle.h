#ifndef GenParticle_h
#define GenParticle_h

#include "TLorentzVector.h"

class GenParticle {
 public:
  enum LeptonType {Muon, Electron} ;
  enum FakeType {notfake, unknown, jet, cjet, bjet, chargemisid};
  enum LooseTight {Loose, Tight, Other};
  
 GenParticle(LeptonType leptonType00, unsigned int leptonIndex00, TLorentzVector& vLepton00, double& eta00, double& chiNdof00, double& D000, double& D0Error00, double& dxy_BS00, double& dz_BS00, int& charge00, FakeType fakeType00, LooseTight looseTight00, double& absTrkIso00, double& relIso00, int& pdgId00, int& mother00)
   : leptonType_(leptonType00), leptonIndex_(leptonIndex00), lorentzVec_(vLepton00), eta_(eta00), chiNdof_(chiNdof00), D0_(D000), D0Error_(D0Error00), dxy_BS_(dxy_BS00), dz_BS_(dz_BS00), charge_(charge00), fakeType_(fakeType00), looseTight_(looseTight00), absTrkIso_(absTrkIso00), relIso_(relIso00), pdgId_(pdgId00), mother_(mother00) {};
  ~GenParticle() {};

  LeptonType leptonType() {return leptonType_; };
  unsigned int ilepton() { return leptonIndex_; };
  
  TLorentzVector& lorentzVec() { return lorentzVec_; }
  
  //void set_relIso(double& reliso) { relIso_ = reliso; }
  void set_relpt_mindr(std::pair<double, double>& thepair) { pair_relpt_mindr_ = thepair; }
  
  int charge() { return charge_; }
  
  double eta() {return eta_; }
  double chiNdof() {return chiNdof_; }
  double D0() { return D0_; }
  double D0Error() { return D0Error_; }
  
  double dxy_BS() { return dxy_BS_; }
  double dz_BS() { return dz_BS_; }
  
  double absTrkIso() { return absTrkIso_; }
  double relIso() { return relIso_; }
  std::pair<double, double>& relpt_mindr() { return pair_relpt_mindr_; }
  
  //void setfakeType(FakeType fakeType) { fakeType_ = fakeType; }
  FakeType fakeType() { return fakeType_; }
  
  LooseTight looseTight() { return looseTight_; }
  void setQuality(LooseTight lt) { looseTight_ = lt; }
  
  int pdfId() {return pdgId_; }
  int mother() {return mother_; }
  //  bool operator<(Lepton& lepton_other) {
  //    return lorentzVec_.Pt() < lepton_other.lorentzVec().Pt();
  //  }

 private:
  LeptonType leptonType_;
  unsigned int leptonIndex_;
  
  TLorentzVector lorentzVec_;
  
  double eta_;
  double chiNdof_;
  double D0_;
  double D0Error_;
  
  double dxy_BS_;
  double dz_BS_;
  
  int charge_;
    
  FakeType fakeType_;
  LooseTight looseTight_;

  double absTrkIso_;
  double relIso_;
  int pdgId_;
  int mother_;
  std::pair<double, double> pair_relpt_mindr_;
};

bool GenPTSorter(GenParticle lep1, GenParticle lep2) {
  return lep1.lorentzVec().Pt() > lep2.lorentzVec().Pt();
}

bool GenIsoSorter(GenParticle lep1, GenParticle lep2) {
  return lep1.relIso() < lep2.relIso();
}

bool GenAbsTrkIsoSorter(GenParticle lep1, GenParticle lep2) {
  return lep1.absTrkIso() < lep2.absTrkIso();
}

//bool LeptonIsolationSorter(Lepton lep1, Lepton lep2) {
//  return (lep1.relIso() < 0.1) > (lep2.relIso() < 0.1);
// true if lep1 earlier-sorted lepton is isolated
//}

#endif
