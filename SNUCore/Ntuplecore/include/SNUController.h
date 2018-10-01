#ifndef SNUController_h
#define SNUController_h

/// standard includes
#include <set>
#include <vector>

// Local include(s):
#include "SNULogger.h"
#include "SNUError.h"
#include "SKTreeFiller.h"

//Forward declaration
class TH1F;
class TChain;

class SNUController  {

 public:
  
  enum _snuversion {none=0,
		    v941=1};

  
  //// constructors
  SNUController();
  virtual ~SNUController();

  enum dataType {NOTSET, data, mc};

  /// Initialise the analysis from the configuration file
  void Initialize() throw( SNUError );  
  void ExecuteCycle()throw( SNUError );

  /// Useful for checking performance
  void GetMemoryConsumption(TString label);
  double GetVirtualMemoryConsumption();
  double GetPhysicalMemoryConsumption();
  void FillMemoryHists(std::string binname);

  int VersionStamp();

  /// List of functions to configure job
  void SetTreeName(TString treename);
  void SetCycleName(TString cyclename);
  void SetSkimName(TString skimname);
  void AddLibraries(TString libname);  
  void SetJobName(TString name);
  void SetTagName(TString name);
  void SetInputList(TString list) throw( SNUError );
  void SetFullInputList(TString list);
  void SetDataType(TString settype);
  void SetLogLevel(TString level);
  void SetName(TString name, Int_t version, TString dir);
  void SetTargetLuminosity(float lumi);
  void SetEffectiveLuminosity(float lumi);
  void SetMCCrossSection(float xsec);
  void SetTotalMCEvents(int total_mc_ev);
  void SetNEventsToProcess(int nevents);
  void SkipEvents(int ev_to_skip);
  void SetOutPutStep( int step);
  void SetDataPeriod(TString period);
  void SetChannel(TString channel);
  void SetInputChain(TChain* ch);
  void SetSNUInput(bool lq);
  void SetUserFlag(TString flag);

  std::string SetNTSNUVersion();

  _snuversion GetSNUVersion(std::string filepath)  throw( SNUError ); 
  bool CheckBranch(std::string ntuple_version, std::string env_version);



  /// Other class functions
  void RunEvent(Long64_t ev);
  void RunNtupleEvent(Long64_t ev);
  void RunTauMode(TString np);
  void RunNonPrompt(TString np);
  void RunChargeFlip(TString cf);

  std::pair< Double_t, Double_t> GetTotalEvents() throw(SNUError);
  float CalculateWeight() throw (SNUError);
  
 private:
  TChain*  chain;
  dataType inputType;
  TString outputLevelString;
  TString CycleName;
  TString skimName;
  TString jobName;
  TString tagName;
  TString treeName;
  TString filelist;
  TString fullfilelist;
  TString completename;
  bool runnp;
  bool runcf;
  bool runtau;
  mutable SNULogger m_logger;
  
  float target_luminosity;
  float sample_crosssection;
  float effective_luminosity;
  double n_total_event;
  int nevents_to_process;
  bool m_isInitialized;
  int n_ev_to_skip;
  
  std::vector<TString> v_libnames;
  std::vector<TString> v_user_flags;
  std::vector<Long64_t> list_to_run;
  Long64_t single_ev;
  Bool_t run_single_event;
  Double_t total_events_beforeskim;
  Double_t total_events_afterskim;
  int output_step;
  TString channel;
  TString k_period;
  bool kSNUInput;
  _snuversion snuversion_lq;

  TH1F* h_timing_hist;
  TH1F* h_virtmemory_hist;
  TH1F* h_physicalmemory_hist;

};
#endif
