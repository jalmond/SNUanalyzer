#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include <iomanip>
#include <sstream>
#include "TSystem.h"

#include <map>

void MakeDataInputFile(TString version=""){
  

  TString def_version = TString(getenv("SNUVERSION"));
  if(!version.Contains("v8") ) version = def_version;

  bool cluster = false;
  TString analysisdir = TString(getenv("HOSTNAME"));
  if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;
  ofstream lumi_file;
  string lfile =  "datasets_snu_SNU_data_" + string(version.Data()) + ".txt";
  if(cluster) lfile =  "datasets_snu_cluster_SNU_data_" + string(version.Data()) + ".txt";

  lumi_file.open(lfile.c_str());
  lumi_file.setf(ios::fixed,ios::floatfield);
  lumi_file.precision(1);
  
  lumi_file << "########################" << endl;
  lumi_file << "### SNUTUPLES ##########" << endl;
  lumi_file << "########################" << endl;
  
  lumi_file << "\n" << endl;
  lumi_file << "#######################################################################################################################################" << endl;
  lumi_file << "# sample              Effective lumi    Path ##########################################################################################" << endl;
  lumi_file << "#######################################################################################################################################" << endl;
  
  
  vector<TString> samples;
  samples.push_back("SingleMuon");
  samples.push_back("SingleElectron");
  samples.push_back("DoubleMuon");
  samples.push_back("DoubleEG");
  samples.push_back("MuonEG");
  samples.push_back("SinglePhoton");
  vector<TString> samples_space;
  samples_space.push_back("    ");
  samples_space.push_back("");
  samples_space.push_back("    ");
  samples_space.push_back("      ");
  samples_space.push_back("        ");
  samples_space.push_back("  ");



  vector<TString> periods;
  periods.push_back("B");
  periods.push_back("C");
  periods.push_back("D");
  periods.push_back("E");
  periods.push_back("F");
  
  TString output="/data8/DATA/SKFlat/v9-4-4/"+version+"/";


  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << samples.at(i)  << samples_space.at(i)<< "            " << periods.at(j) <<   "        " <<  output +  samples.at(i) << "/period" <<  periods.at(j) << "/" << endl;
    }
    //lumi_file << ""<< endl;
  }

  
  TString SKTreeOutput="/data8/DATA/SnuNtuples/"+ version + "/SKTrees/Data/";
  if(cluster) SKTreeOutput="/data4/LocalNtuples/SKTrees13TeV/"+ version + "/SKTrees/Data/";

  lumi_file << ""<< endl;
  lumi_file << "#### Single lepton skims: SKTrees" << endl;
  lumi_file << ""<< endl;
  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << "SK" << samples.at(i)  << samples_space.at(i) << "            " << periods.at(j) <<   "        " <<  SKTreeOutput +  samples.at(i) << "/period" <<  periods.at(j) << "/" << endl;
    }
    //lumi_file << ""<< endl;
  }

  lumi_file << ""<< endl;
  lumi_file << "#### Single dilepton skims: SKTrees" << endl;
  lumi_file << ""<< endl;
  

  TString SKTreeOutputDiLep="/data7/DATA/SnuNtuples/"+ version + "/SKTrees/DataDiLep/";
  if(cluster) SKTreeOutputDiLep="/data4/LocalNtuples/SKTrees13TeV/"+ version + "/SKTrees/DataDiLep/";

  
  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << "SK" << samples.at(i) << "_dilep"  << samples_space.at(i) <<"            " << periods.at(j) <<   "        " <<  SKTreeOutputDiLep +  samples.at(i) << "/period" <<  periods.at(j)<< "/" << endl;
    }
    //lumi_file << ""<< endl;
  }
  
  TString SKTreeOutputHNDiLep="/data7/DATA/SnuNtuples/"+ version + "/SKTrees/DataHNDiLep/";


  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << "SK" << samples.at(i) << "_hndilep"  << samples_space.at(i) <<"            " << periods.at(j) <<   "        " <<  SKTreeOutputHNDiLep +  samples.at(i) << "/period" <<  periods.at(j)<< "/" << endl;
    }
    //lumi_file << ""<< endl;                                                                                                                                                     
  }

  
  TString SKTreeOutputHNFake="/data7/DATA/SnuNtuples/"+ version + "/SKTrees/DataHNFake/";
  
  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << "SK" << samples.at(i) << "_hnfake"  << samples_space.at(i) <<"            " << periods.at(j) <<   "        " <<  SKTreeOutputHNFake +  samples.at(i) << "/period" <<  periods.at(j)<< "/" << endl;
    }
    //lumi_file << ""<< endl;                                                                                                                                                   \
                                                                                                                                                                                 
  }

 


  TString SKTreeOutputTriLep="/data7/DATA/SnuNtuples/"+ version + "/SKTrees/DataTriLep/";
  if(cluster) SKTreeOutputTriLep="/data4/LocalNtuples/SKTrees13TeV/"+ version + "/SKTrees/DataTriLep/";

  for(unsigned int i = 0 ; i < samples.size() ; i++){
    for(unsigned int j = 0 ; j < periods.size() ; j++){
      lumi_file << "SK" << samples.at(i) << "_trilep"  << samples_space.at(i) <<"            " << periods.at(j) <<   "        " <<  SKTreeOutputTriLep +  samples.at(i) << "/period" <<  periods.at(j)<< "/" << endl;
    }
    //lumi_file << ""<< endl;
  }


  string lqdir = getenv("LQANALYZER_DIR");
  TString user = TString(getenv("USER"));

  string lfile2 =   lqdir+ "/LQRun/txt/datasets_snu_SNU_data_" + string(version.Data()) + ".txt";
  if(cluster) lfile2 =   lqdir+ "/LQRun/txt/Cluster/datasets_snu_cluster_SNU_data_" + string(version.Data()) + ".txt";

  if(user.Contains("jalmond")){
    if(!cluster)gSystem->Exec(("cp " + lfile + "  " + getenv("LQANALYZER_DATASETFILE_DIR")).c_str());     
    else gSystem->Exec(("cp " + lfile + "  " + getenv("LQANALYZER_DATASETFILE_DIR")).c_str());

  }
  gSystem->Exec(("mv " + lfile +" " + lfile2).c_str());

    
  return;
  
}
