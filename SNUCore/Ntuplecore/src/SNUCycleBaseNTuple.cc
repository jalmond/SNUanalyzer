/// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>

/// STL includes
#include <iostream>
#include <iomanip>

//local  includes
#include "SNUCycleBaseNTuple.h"


ClassImp( SNUCycleBaseNTuple);

SNUCycleBaseNTuple::SNUCycleBaseNTuple() : SNUCycleBaseBase(), m_outputFile(0),m_outputTrees(),m_outputVarPointers(), k_isdata(false) ,k_running_nonprompt(false), k_running_chargeflip(false),k_running_taudecays(false), k_sample_name(""), k_tag_name(""),k_classname(""),k_skim(""), sample_entries(-999), output_interval(10000), events_to_process(-1), k_mcperiod(0) {


 
}


SNUCycleBaseNTuple::~SNUCycleBaseNTuple(){
  DeleteInputVariables();
}



void SNUCycleBaseNTuple::DeleteInputVariables() {

  for( std::list< TObject* >::iterator it = m_inputVarPointers.begin();
       it != m_inputVarPointers.end(); ++it ) {
    delete ( *it );
  }
  m_inputVarPointers.clear();

  return;
}

void SNUCycleBaseNTuple::CreateOutputTrees(TFile* outputFile, TString name, TString title){
  
  m_logger << INFO  << "Creating Output Trees" << SNULogger::endmsg;

  // Clear the vector of output trees:
  m_outputTrees.clear();
  
  // Clear the vector of output variable pointers:                             
  m_outputVarPointers.clear();

  // Open output file / create output trees 
  gROOT->cd();

  // Create all the regular output trees, but don't create any branches in them
  // just yet.
  //
  
  const Int_t branchStyle = 1;
  const Int_t autoSave = 10000000;
  
  TTree* tree = new TTree(name,title);
  tree->SetAutoSave( autoSave );
  TTree::SetBranchStyle( branchStyle );
  
  m_outputTrees.push_back( tree );
  
  if( outputFile ) {
    tree->SetDirectory( outputFile );
  }
  
  return;
}


void SNUCycleBaseNTuple::SetOutPutStep(int step){

  if(step != output_interval) m_logger << INFO << "Changing default value for output interval from every 10000 events to every " << step << " events" << SNULogger::endmsg;
  output_interval = step;  
}

void SNUCycleBaseNTuple::SetNEventsToProcess(int nentries){
  events_to_process = nentries;
}


void SNUCycleBaseNTuple::SetDataType( bool type){  
  k_isdata = type;
}

void SNUCycleBaseNTuple::SetMCPeriod( int period){
  m_logger << INFO << "SetMCPeriod : " << k_mcperiod << SNULogger::endmsg;
  k_mcperiod= period;
}

void SNUCycleBaseNTuple::SetTauStatus( bool type){
  k_running_taudecays = type;
}


void SNUCycleBaseNTuple::SetNPStatus( bool type){
  k_running_nonprompt = type;
}

void SNUCycleBaseNTuple::SetCFStatus( bool type){
  k_running_chargeflip = type;
}

void SNUCycleBaseNTuple::SetSampleName( TString name){
  k_sample_name = name;
}
void SNUCycleBaseNTuple::SetTagName( TString name){
  k_tag_name= name;
}
void SNUCycleBaseNTuple::SetAnalyzerClassName( TString c){
  k_classname= c;
}

void SNUCycleBaseNTuple::SetDataChannel( TString channel){
  k_channel= channel;
}
void SNUCycleBaseNTuple::SetSkimName( TString skim){
  k_skim= skim;
}



void SNUCycleBaseNTuple::SetNSampleEvents(double nev){
  sample_entries = nev;
}

TFile* SNUCycleBaseNTuple::GetOutputFile(){
  return m_outputFile;
}
  
void SNUCycleBaseNTuple::MakeOutPutFile(TString outfile, TString treename){
  
  if(!m_outputFile){
    m_logger << INFO << "Creating " << outfile << SNULogger::endmsg;  
    m_outputFile = TFile::Open(outfile, "RECREATE");
  }else {
    m_logger << WARNING << "Output file created already. Will not create again" << SNULogger::endmsg;
  }
  
  if(!treename.Contains("NOTREE"))  CreateOutputTrees(m_outputFile, treename , "");
  
}

void SNUCycleBaseNTuple::FillOutTree(){

  int nbytes = 0;
  for( std::vector< TTree* >::iterator tree = m_outputTrees.begin();
       tree != m_outputTrees.end(); ++tree ) {
    nbytes = ( *tree )->Fill();
    if( nbytes < 0 ) {
      REPORT_ERROR( "Write error occured in tree \""
		    << ( *tree )->GetName() << "\"" );
    } else if( nbytes == 0 ) {
      m_logger << WARNING << "No data written to tree \""
	       << ( *tree )->GetName() << "\"" << SNULogger::endmsg;
    }
  }
}

void SNUCycleBaseNTuple::WriteCycleHists(TH1F* h_timing, TH1F* hmem1, TH1F* hmem2){

  // Remember which directory we were in: 
  TDirectory* savedir = gDirectory;
  
  TDirectory* cycledir = m_outputFile->mkdir("CycleInfo");
  m_outputFile->cd( cycledir->GetName());
  h_timing->Write();
  hmem1->Write();
  hmem2->Write();

  // Go back to the original directory:                                                                                                                              
  gDirectory = savedir;
  
}

void  SNUCycleBaseNTuple::CloseFiles(){
  
  m_logger << INFO << "Closing output file  " << m_outputFile->GetName() << SNULogger::endmsg;
  m_outputFile->SaveSelf( kTRUE ); /// is this needed 
  m_outputFile->Close();
  delete m_outputFile;
  m_outputFile = 0;
  
}

void SNUCycleBaseNTuple::SaveOutputTrees( TDirectory* /*output*/ ) {

  // Remember which directory we were in:
  TDirectory* savedir = gDirectory;

  // Save each regular output tree:

  for( std::vector< TTree* >::iterator tree = m_outputTrees.begin();
       tree != m_outputTrees.end(); ++tree ) {
    TDirectory* dir = ( *tree )->GetDirectory();
    if( dir ) dir->cd();
    ( *tree )->Write();
    ( *tree )->AutoSave();
    ( *tree )->SetDirectory( 0 );
    delete ( *tree );
  }

  // Go back to the original directory:                                                                    
  gDirectory = savedir;

  return;
}

TTree* SNUCycleBaseNTuple::GetOutputTree( const char* treeName ) const{
  
  // Look for such output tree:                                                
  //                                                                           
  m_logger << INFO << "Getting output Tree " << SNULogger::endmsg;
  TString tname( treeName );
  for( std::vector< TTree* >::const_iterator it = m_outputTrees.begin();
       it != m_outputTrees.end(); ++it ) {
    if( *it ) {
      if( tname == ( *it )->GetName() ) {
	return *it;
      }
    }
  }
  m_logger << INFO << "Asked for output tree when none was initially created. " << SNULogger::endmsg;
  return 0;
}


const char* SNUCycleBaseNTuple::RootType( const char* typeid_type ){

  switch( typeid_type[ 0 ] ) {

  case 'c':
    return "B";
    break;
  case 'h':
    return "b";
    break;
  case 's':
    return "S";
    break;
  case 't':
    return "s";
    break;
  case 'i':
    return "I";
    break;
  case 'j':
    return "i";
    break;
  case 'f':
    return "F";
    break;
  case 'd':
    return "D";
    break;
  case 'x':
    return "L";
    break;
  case 'y':
    return "l";
    break;
  case 'b':
    return "O";
    break;

  }

  return "";
}
