#ifndef FakeRateCalculator_FinalMu_h
#define FakeRateCalculator_FinalMu_h

#include "AnalyzerCore.h"


class FakeRateCalculator_FinalMu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_FinalMu();
  ~FakeRateCalculator_FinalMu();

  /// Functions from core
  virtual void BeginCycle() throw( SNUError );
  virtual void BeginEvent()throw( SNUError );
  virtual void ExecuteEvents()throw( SNUError );
  virtual void EndCycle()throw( SNUError );
  virtual void ClearOutputVectors()throw( SNUError );
  
  void InitialiseAnalysis() throw( SNUError );
  void MakeHistograms();


  void GetFakeRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el,TString tightlabel,  std::vector<snu::KJet> jets,std::vector<snu::KJet> alljets, TString tag, double w, bool makebasicplots);
  void GetHSTRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el, std::vector<snu::KJet> jets, TString tag);
  
  void MakeMCFakeratePlots(TString label, bool pass_single_trigger, std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, double w);
  void MakePlotsMCAwaJetPt(TString label,float awayjetptcut, std::vector<snu::KElectron>, std::vector<snu::KJet> jets , std::vector<snu::KJet> alljets, double w);

  
  void MakeMCPlots(TString label, snu::KElectron electron, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, double w);
  float GetPrescale( std::vector<snu::KElectron> electrons, bool passlowest, bool passlow, bool passhigh, bool highest, bool passshighest, float fake_total_lum);
  bool  UseEvent(std::vector<snu::KElectron> electrons,  vector<snu::KJet> jets, float pcut,float precale_weight, float w);
  bool  UseEventAllMT(std::vector<snu::KElectron> electrons,  vector<snu::KJet> jets, float pcut,float precale_weight, float w);



  void MakeFakeRatePlots(TString label, TString eltag,   std::vector<snu::KElectron> tightelectrons,std::vector<snu::KElectron> electrons,   std::vector <snu::KJet> jets, std::vector<snu::KJet> alljets,  float precale_w, float w, bool makebasicplots);
  void MakeDXYFakeRatePlots(TString label, TString eltag,  std::vector<snu::KElectron> electrons,   std::vector <snu::KJet> jets, std::vector<snu::KJet> alljets,  float precale_w, float prescale_diel, float w);



  void MakeSingleElectronCRPlots(TString looseid, TString eltag, TString tightid, float w, bool usepujetid);

  void GetFakeRateAndPromptRates(std::vector<snu::KElectron> electrons, TString eltag, std::vector<snu::KElectron> tightelectrons, float w, bool usepujetid, bool runall);

  void GetFakeRateAndPromptRatesPerPeriod(TString looseid, TString eltag, TString tightid, float w, bool usepujetid, bool runall);
private:
  
  
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;

  int n_17_pass;
  int n_17_17_jet_pass;
  int n_17_jet_pass;
  
  ClassDef ( FakeRateCalculator_FinalMu, 1);
};
#endif
