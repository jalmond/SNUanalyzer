// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: SNUExampleAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: SNUCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in SNUCore/core classes
ClassImp (ExampleAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ExampleAnalyzer::ExampleAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ExampleAnalyzer");
  
  Message("In ExampleAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ExampleAnalyzer::InitialiseAnalysis() throw( SNUError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << SNULogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:


  return;
}


void ExampleAnalyzer::ExecuteEvents()throw( SNUError ){

  /// Apply the gen weight 
  if(!IsData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << SNULogger::endmsg;
  m_logger << DEBUG << "IsData = " << IsData << SNULogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  if(IsData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);


   if(!PassMETFilter()) return;     /// Initial event cuts : 
   FillCutFlow("EventCut", weight);

   /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
   
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              

   float pileup_reweight=(1.0);
   if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
     
   
   TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
   TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
   TString dimuon_trigmuon_trig3="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
   TString dimuon_trigmuon_trig4="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
   // Now you should do an OR of 4 triggers 
 
   vector<TString> trignames;
   trignames.push_back(dimuon_trigmuon_trig1);
   trignames.push_back(dimuon_trigmuon_trig2);
   trignames.push_back(dimuon_trigmuon_trig3);
   trignames.push_back(dimuon_trigmuon_trig4);


   std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_NOCUT");
   std::vector<snu::KElectron> CFelectrons = GetElectrons(false,false, "ELECTRON_HN_TIGHTv4");
   /*
     
   std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_NOCUT);  ... WONT WORK
   std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_NOCUT");               ... WILL WORK  
   
   std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);  ... WILL WORK  
   std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_POG_TIGHT");                ... WILL WORK  
   
   */

   //   std::vector<snu::KElectron> electrons2 =  GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0);

   std::vector<snu::KJet> jets =   GetJets("JET_HN");
   int nbjet = NBJet(GetJets("JET_HN"));
   std::vector<snu::KMuon> muons =GetMuons("MUON_HN_TIGHT",false); 

   if(muons.size() > 0) cout << "muon reliso = " << muons[0].RelIso04() << endl;
   bool trig_pass= PassTriggerOR(trignames);


   mcdata_correction->CorrectMuonMomentum(muons,eventbase->GetTruth()); /// CorrectMuonMomentum(muons);  will also work as Funcion in AnalyzerCore just calls mcdata_correction function
   
   double ev_weight = weight;
   if(!IsData){
     //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
   }

   if(jets.size() > 3){
     if(nbjet > 0){
       if(muons.size() ==2) {
	 if(electrons.size() >= 1){
	   cout << "electrons is tight = " << electrons.at(0).PassTight() << endl;
	   if(!SameCharge(muons)){
	     if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 10.){
	       if(eventbase->GetEvent().PFMET() > 30){
		 if(trig_pass){
		   FillHist("Massmumu", GetDiLepMass(muons), ev_weight, 0., 200.,400);
		   FillHist("Massmumu_zoomed", GetDiLepMass(muons), ev_weight, 0.,50.,200);
		   FillCLHist(sighist_mm, "DiMuon", eventbase->GetEvent(), muons,electrons,jets, ev_weight);
		 }
	       }
	     }
	   }
	 }
       }
     }
   }

   	    
   float cf_weight = -999.;
   if(CFelectrons.size() == 2){
     if(CFelectrons.at(0).Charge() != CFelectrons.at(1).Charge()){//CF estimation is from OS dielectrons
       cf_weight = GetCFweight(0,CFelectrons, true, "ELECTRON_HN_TIGHTv4");//put in syst=1,0,-1,  electronColl vector, apply sf, electron ID
         //scale factor up downs with syst = 1, -1

       std::vector<snu::KElectron> CFelectrons_shifted = ShiftElectronEnergy(CFelectrons, "ELECTRON_HN_TIGHTv4", true);//apply pt shift considering brem radiation after getting CF weights
       if(CFelectrons_shifted.at(1).Pt()>20){//apply pt cuts again after shifting pt
         FillHist("chargeflipped_Z_mass", (CFelectrons.at(0)+CFelectrons.at(1)).M(),cf_weight, 50., 130., 80);
         FillHist("chargeflipped_leading_lepton", CFelectrons.at(0).Pt(), cf_weight, 0., 200., 200);
   
       }
     }
   }

   
   return;
}// End of execute event loop
  


void ExampleAnalyzer::EndCycle()throw( SNUError ){
  
  Message("In EndCycle" , INFO);

}


void ExampleAnalyzer::BeginCycle() throw( SNUError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: SNUTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "SNUTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

ExampleAnalyzer::~ExampleAnalyzer() {
  
  Message("In ExampleAnalyzer Destructor" , INFO);
  
}


void ExampleAnalyzer::BeginEvent( )throw( SNUError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ExampleAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ExampleAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ExampleAnalyzer::ClearOutputVectors() throw(SNUError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



