#!/bin/sh
### sets all configurable variables to defaul values

source    /data1/SNUAnalyzer_rootfiles_for_analysis/SnutupleConfig/${SNUVERSION}.sh
cp $ANALYZER_DATASETFILE_DIR/datasets_snu* $ANALYZER_DIR/SNURun/txt/
cp $ANALYZER_DATASETFILE_DIR/list_all_mc*  $ANALYZER_DIR/SNURun/txt/


declare -a list_of_skims=("SKFLAT" "SKTree_NoSkim" "SKTree_LeptonSkim" "SKTree_DiLepSkim" "SKTree_TriLepSkim" "SKTree_SSLepSkim" "NoCut" "Lepton" "DiLep")
declare -a list_of_sampletags=("ALL" "DATA" "MC" "DoubleEG" "DoubleMuon" "MuonEG" "SingleMuon" "SinglePhoton" "SingleElectron" "SingleLepton" "DoubleMuon_CF")

#declare -a queueoptions=("allq" "fastq" "longq" "node1" "node2" "node3" "node4" "node5" "node6" "None")
declare -a queueoptions=("allq" "fastq" "longq"  "None")  

python $ANALYZER_DIR/python/CheckSelection.py
if [[ ! -d  /data2/SNU_SKTreeOutput/${USER}/GoodSelection/ ]]; then
    echo "Fix selection file and rerun"
    exit 1
fi

##### New for sktreemaker only
logger=""

###### SET WHAT JOBS TO RUN
runMC=false
runDoubleMuon=false
runDoubleElectron=false
runSingleMuon=false
runElectronMuon=false
runSingleElectron=false
runSinglePhoton=false

## RUN PARAMETERS

job_data_lumi="ALL"   ### ALL = "BtoG"  ###  "C" = period C only   "ALL"  = period C+D
job_submitallfiles="false"
job_logstep=1000
job_loglevel="INFO"
job_njobs=2
job_skim="SKTree_LeptonSkim"
changed_skim=false
job_output_dir=""
job_tmp_filename="None"

RUNSIG="false"
make_sktrees="False"
GetOutPutDir="False"
run_validation="False"
idname=""
queuename="None"
object=""
idname2=""
job_run_bkg="False"
submit_draw="False"
submit_sk_message=""
submit_skflag=""
submit_skinput=true
submit_dir_tag="NULL"
changed_skinput=false
changed_submit_version_tag=false
changed_job_output_dir=false
changed_job_loglevel=false
changed_job_data_lumi=false
changed_job_njobs=false 
submit_file_tag=""
submit_file_list=""
submit_version_tag=${SNUVERSION}
submit_sampletag=""
submit_snuvlist=""
submit_searchlist=""
submit_analyzer_name=""
set_submit_analyzer_name=false
request_sample=""
submit_skim=""
submit_snu_tag=""
submit_snu_tag2=""

######## NEW FOR TAG v7-6-3.2
job_nevents=-1
job_nskip=-1
job_run_taus=False
job_run_fake=False
job_run_flip=False
check_all_snuversions=false
set_submit_file_tag=false
set_submit_file_list=false
set_submit_sampletag=false
set_sktreemaker_debug=false



FLATSNU_MC="/data2/DATA/snutoflat/MC/"
SKTREE_MC="/data2/SnuNtuples/"
if [ $HOSTNAME == "cmscluster.snu.ac.kr" ];
then
    TXTPATH=${ANALYZER_RUN_PATH}"/txt/Cluster/datasets_snu_cluster_"
    FLATSNU_MC="/data4/DATA/FlatSnuuples/MC/"
    SKTREE_MC="/data4/LocalNtuples/SKTrees13TeV/"
fi

is_mc=""

### Get predefined lists
source ${ANALYZER_DIR}/SNURun/txt/list_all_mc_${submit_version_tag}.sh
### setup list of samples and other useful functions
TXTPATH=${ANALYZER_RUN_PATH}"/txt/datasets_snu_nonsig_"
source submit_setup.sh

OLDTXTPATH=${ANALYZER_RUN_PATH}"/txt/datasets_snu_"


if [[ $submit_sampletag  == "" ]];
then
    is_mc="true"
fi



linkdef_filepath=/SNUAnalysis/Analyzers/include/SNUAnalysis_LinkDef.h

if [[ $make_sktrees == "True" ]];
    then
    linkdef_filepath=/SNUAnalysis/SKTreeMaker/include/SKTreeMaker_LinkDef.h
fi

if [[ $run_validation == "True" ]];
    then
    linkdef_filepath=/SNUAnalysis/Validation/include/Validation_LinkDef.h
fi

############## Check flags for fake/flip analysis
if [[ $job_run_fake != "False" ]];
    then
    if [[ $job_run_fake != "True" ]];
        then
	echo "SNUanalyzer::sktree :: ERROR :: Wrong setting for -fake. use -fake 'True'. It is false by default".
        exit 1
  fi
fi
if [[ $job_run_taus != "False" ]];
    then
    if [[ $job_run_taus != "True" ]];
        then
        echo "SNUanalyzer::sktree :: ERROR :: Wrong setting for -fake. use -fake 'True'. It is false by default".
        exit 1
  fi
fi

if [[ $job_run_flip != "False" ]];
    then
    if [[ $job_run_flip != "True" ]];
        then
        echo "SNUanalyzer::sktree :: ERROR :: Wrong setting for -flip. use -flip 'True'. It is false by default".
        exit 1
  fi
fi


########### safe guard wrong snuversion setting
if [[ $check_all_snuversions == "true" ]];
    then
    submit_version_tag=""
    if [[ $changed_submit_version_tag == "true" ]];
        then
        echo "SNUanalyzer::sktree :: ERROR :: Snuversion was set with -c, while you chose to look at all snuversions wit -ac option. This is not possible."
        exit 1
     fi
fi





if [[ $set_submit_analyzer_name == "true" ]];
    then
    if [[ $submit_analyzer_name ==  "" ]];
	then
	echo "No analyzer set: set with -a"
	echo "List of Classes are:"
	while read line
	  do
	  if [[ $line == *"C++ class"* ]];
	      then
	      if [[ $line != *"AnalyzerCore"* ]];
		  then
		  sline=$(echo $line | head -n1 | awk '{print $5}')
		  suffix="+;"
		  sline=${sline%$suffix}
		  echo $sline
	      fi
	  fi
	done < ${ANALYZER_DIR}/$linkdef_filepath
	exit 1
    fi
fi


if [[ $submit_analyzer_name !=  "" ]];
then
    analyzer_found=false
    while read line
    do
	if [[ $line == *"C++ class"* ]];
        then
	    if [[ $line != *"AnalyzerCore"* ]];
	    then
		sline=$(echo $line | head -n1 | awk '{print $5}')
		suffix="+;"
		sline=${sline%$suffix}
		if [[ $sline == $submit_analyzer_name ]];
		then
		    analyzer_found=true
		fi
	    fi
	fi
    done < ${ANALYZER_DIR}/$linkdef_filepath
    
    if [[ $analyzer_found == "false" ]];
    then 
	echo ""
	echo "SNUanalyzer::sktree :: ERROR :: Analyzer is invalid: 'sktree -a <analyzer/class>'"
	echo "Allowed Analyzers are:"
	while read line
	do
	    if [[ $line == *"C++ class"* ]];
	    then
		if [[ $line != *"AnalyzerCore"* ]];
		then
		    sline=$(echo $line | head -n1 | awk '{print $5}')
		    suffix="+;"
		    sline=${sline%$suffix}
			echo $sline
		fi
	    fi
	done < ${ANALYZER_DIR}/$linkdef_filepath
	exit 1
    fi
fi

function mergeoutput
{
    if [[ $job_cycle != *"SKTreeMaker"* ]];
        then


	if [[ $job_data_lumi  == "ALL" ]];
            then
	    output_file_skim_tag=$1
	    
	    if [[ $job_skim == "FLATSNU" ]];
		then
		output_file_skim_tag=$output_file_skim_tag"_snu_"$submit_version_tag

	    fi
	    if [ $job_skim == "SKTree_NoSkim" ] || [ $job_skim == "NoCut" ];
		  then
                output_file_skim_tag="SK"$output_file_skim_tag"_nocut_snu_"$submit_version_tag
            fi
	    
	    if [ $job_skim == "SKTree_LeptonSkim" ] || [ $job_skim == "Lepton" ];
		then
		output_file_skim_tag="SK"$output_file_skim_tag"_snu_"$submit_version_tag
	    fi
	    
	    if [ $job_skim == "SKTree_DiLepSkim" ] || [ $job_skim == "DiLep" ];
		then
		output_file_skim_tag="SK"$output_file_skim_tag"_dilep_snu_"$submit_version_tag
	    fi

	    if [[ $job_skim == "SKTree_TriLepSkim" ]] ;
		then
		output_file_skim_tag="SK"$output_file_skim_tag"_trilep_snu_"$submit_version_tag
            fi
	    if [[ $job_skim == "SKTree_SSLepSkim" ]] ;
	    then
                output_file_skim_tag="SK"$output_file_skim_tag"_sslep_snu_"$submit_version_tag
            fi
	    
	    echo "############################################################################################################################################################"

	    echo "MERGED DATA: "
	    echo "Command:"
            #echo "source hadd.sh "${outputdir_data}" "${job_cycle}"_data_snu_"${submit_version_tag}".root "${outputdir_data}${job_cycle}"'*'"$output_file_skim_tag"'*'"
            #echo "############################################################################################################################################################"
	    #echo ""
            #source hadd.sh ${outputdir_data} ${job_cycle}_data_snu_${submit_version_tag}.root ${outputdir_data}${job_cycle}'*'${output_file_skim_tag}'*'
	    echo ""
            echo "merged output sent to -----> "${outputdir_data}${job_cycle}"_data_snu_"${submit_version_tag}".root "
	    
            #mv  ${outputdir_data}/${job_cycle}_data_snu_${submit_version_tag}.root  ${outputdir_mc}/${job_cycle}_data_$1_snu_${submit_version_tag}.root

            echo ""
        fi
    fi

}

function mergefake
{
            output_file_skim_tag=$1
	    outname="nonprompt"

            if [[ $job_skim == "FLATSNU" ]];
                then
                output_file_skim_tag=$output_file_skim_tag"_snu_"$submit_version_tag
		outname=$outname"_snu_"$submit_version_tag
            fi
            if [ $job_skim == "SKTree_NoSkim" ] || [ $job_skim == "NoCut" ];
                  then
                output_file_skim_tag="SK"$output_file_skim_tag"_nocut_snu_"$submit_version_tag
		outname="SK"$outname"_nocut_snu_"$submit_version_tag
            fi

            if [ $job_skim == "SKTree_LeptonSkim" ] || [ $job_skim == "Lepton" ];
                then
                output_file_skim_tag="SK"$output_file_skim_tag"_snu_"$submit_version_tag
		outname="SK"$outname"_snu_"$submit_version_tag

            fi

            if [ $job_skim == "SKTree_DiLepSkim" ] || [ $job_skim == "DiLep" ];
                then
                output_file_skim_tag="SK"$output_file_skim_tag"_dilep_snu_"$submit_version_tag
		outname="SK"$outname"_dilep_snu_"$submit_version_tag

            fi

            if [[ $job_skim == "SKTree_TriLepSkim" ]] ;
                then
                output_file_skim_tag="SK"$output_file_skim_tag"_trilep_snu_"$submit_version_tag
		outname="SK"$outname"_trilep_snu_"$submit_version_tag

            fi
	    if [[ $job_skim == "SKTree_SSLepSkim" ]] ;
                then
                output_file_skim_tag="SK"$output_file_skim_tag"_sslep_snu_"$submit_version_tag
                outname="SK"$outname"_sslep_snu_"$submit_version_tag

            fi

	    #echo hadd.sh ${outputdir_np} ${job_cycle}_${output_file_skim_tag}.root ${outputdir_np}${job_cycle}'*'${output_file_skim_tag}'*'
	    if [[ $job_data_lumi  == "ALL" ]];
	    then
		#source hadd.sh ${outputdir_np} ${job_cycle}_${output_file_skim_tag}.root ${outputdir_np}${job_cycle}'*'${output_file_skim_tag}'*'
		#mv  ${outputdir_np}/${job_cycle}_${output_file_skim_tag}.root  ${outputdir_mc}/${job_cycle}_$1_${outname}.root
		echo "Merged data will be sent to ---> " +  ${outputdir_mc}/${job_cycle}_$1_${outname}.root
	    fi
	    if [[ $job_data_lumi  == $snudatatag ]];
	    then
                #source hadd.sh ${outputdir_np} ${job_cycle}_${output_file_skim_tag}.root ${outputdir_np}${job_cycle}'*'${output_file_skim_tag}'*'
		#mv  ${outputdir_np}/${job_cycle}_${output_file_skim_tag}.root  ${outputdir_mc}/${job_cycle}_$1_${outname}.root
		echo "Merged data will be sent to ---> " +  ${outputdir_mc}/${job_cycle}_$1_${outname}.root
            fi
}

if [[ $changed_skim == "true" ]];
then
    skimisok=false
    for iskim in  ${list_of_skims[@]};
    do
	if [[ $job_skim == $iskim ]];
	then
	    skimisok=true
	fi
    done
    if [[ $skimisok == "false" ]];
    then
	echo "SNUanalyzer::sktree :: ERROR :: Invalid skim -s <skim> "
	echo "Allowed values are:"
	  for iskim in  ${list_of_skims[@]};
	  do
	      echo $iskim
	  done
	  exit 1
    fi
fi

if [[ $submit_file_tag  != ""  ]];
    then
    ARG_SINGLE_FILE="FULLLISTOFSAMPLES"
    if [[ $job_skim == "SKTree_NoSkim" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESNOCUT"
    fi
    if [[ $job_skim == "NoCut" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESNOCUT"
    fi
    if [[ $job_skim == "SKTree_LeptonSkim" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESLEPTON"
    fi
    if [[ $job_skim == "Lepton" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESLEPTON"
    fi
    if [[ $job_skim == "SKTree_DiLepSkim" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESDILEP"
    fi


    if [[ $job_skim == "SKTree_TriLepSkim" ]];then
        ARG_SINGLE_FILE="FULLLISTOFSAMPLESTRILEP"
    fi
    if [[ $job_skim == "SKTree_SSLepSkim" ]];then
        ARG_SINGLE_FILE="FULLLISTOFSAMPLESSS"
    fi

    if [[ $job_skim == "DiLep" ]];then
	ARG_SINGLE_FILE="FULLLISTOFSAMPLESDILEP"
    fi

    eval test_single_file=(\${$ARG_SINGLE_FILE[@]})
    file_tag_exists=false
    for all_files in ${test_single_file[@]};
      do
	if [[ $all_files == $submit_file_tag ]];
	then
	    file_tag_exists=true
	fi
    done
    
    if [[ $submit_sampletag != "MC" ]];
	then
	if [[ $submit_file_list  != ""  ]];
	    then
	    echo "SNUanalyzer::sktree :: ERROR :: Running with 'sktree -i and -list is not valid'"
	    exit 1
	fi
    fi
    list_isdata=true
    if [[ $submit_sampletag == "MC" ]];
	then
	list_isdata=false
    fi
    if [[ $list_isdata == "false" ]];
	then
	
	if [[ $submit_sampletag != ""  ]];
	    then
	     echo "SNUanalyzer::sktree :: ERROR :: Running with 'sktree -i <samplename> and -S <samplelist> is not valid unless samplelist is data'"
	    exit 1
	fi
    fi


fi

if [[ $submit_file_tag  == ""  ]];
    then
    if [[ $submit_file_list  == ""  ]];
	then
	if [[ $submit_sampletag == ""  ]];
            then
	    
	    if [[ $set_submit_sampletag  == "true"  ]];
		then
		echo "SNUanalyzer::sktree :: ERROR :: no input of sktree -s <SAMPLE>"
		echo  "Options are:"
		for isample in ${list_of_sampletags[@]};
		  do
		  echo $isample
		done
		exit 1
	    fi

	    if [[ $set_submit_file_tag  == "true"  ]];
                then
		  echo "SNUanalyzer::sktree :: ERROR :: no input of sktree -i <samplename> " 
		  submit_skim=$job_skim
		  submit_searchlist=""
		  submit_snuvlist=$submit_version_tag
		  runlist
		  exit 1
	    fi

	    if [[ $set_submit_file_list  == "true"  ]];
                then
		echo "SNUanalyzer::sktree :: ERROR :: no input of sktree -list <list> "
		submit_searchlist=""
		submit_snuvlist=$submit_version_tag
		rungroupedlist
		echo "SNUanalyzer::sktree :: HELP :: to check content of list run 'sktree -g list'"
		echo "SNUanalyzer::sktree :: HELP :: You can add your own list to "$ANALYZER_DIR"/SNURun/txt/list_user_mc.sh"
		exit 1
            fi
	    
            echo "SNUanalyzer::sktree :: ERROR :: No input files were specified"
	    echo "Use either:"
	    echo "sktree -S <SAMPLELIST>   : run 'sktree -h' for list of options"
	    echo "sktree -i <SAMPLENAME>     : run 'sktree -l' or 'sktree -L' for list of options "
	    echo "sktree -list <SAMPLEARRAY> : run 'sktree -g' for list of options  "
            exit 1
        fi

    fi
fi


    


if [[ $submit_file_list  != ""  ]];
    then
    valid_list_samples=0
    ARGTEST=$submit_file_list
    eval test_array=(\${$ARGTEST[@]})
    array_exists=false
    counter=${#test_array[@]}
	if [[ $counter -eq 0 ]];
            then
	    echo "SNUanalyzer::sktree :: ERROR :: Input list is invalid: 'sktree -list list'"
	    echo "SNUanalyzer::sktree :: HELP :: Check lists using 'sktree -g'"
	    exit 1
	fi
	
	for part_of_list in ${test_array[@]}; do
	    if [[ $part_of_list == "" ]];
		then
		echo "SNUanalyzer::sktree :: ERROR :: Check input list. You have missing sample ''."
                exit 1
	    fi
	    ARG_SINGLE_FILE="FULLLISTOFSAMPLES"
	    if [[ $job_skim == "SKTree_NoSkim" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESNOCUT"
	    fi
	    if [[ $job_skim == "NoCut" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESNOCUT"
	    fi
	    if [[ $job_skim == "SKTree_LeptonSkim" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESLEPTON"
	    fi
	    if [[ $job_skim == "Lepton" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESLEPTON"
	    fi
	    if [[ $job_skim == "SKTree_DiLepSkim" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESDILEP"
	    fi

	    if [[ $job_skim == "SKTree_TriLepSkim" ]];then
                ARG_SINGLE_FILE="FULLLISTOFSAMPLESTRILEP"
            fi
	    if [[ $job_skim == "SKTree_SSLepSkim" ]];then
                ARG_SINGLE_FILE="FULLLISTOFSAMPLESSS"
            fi
	    if [[ $job_skim == "DiLep" ]];then
		ARG_SINGLE_FILE="FULLLISTOFSAMPLESDILEP"
	    fi

	    
	    eval test_single_file=(\${$ARG_SINGLE_FILE[@]})
	    file_tag_exists=false
	    for all_files in ${test_single_file[@]};
	      do
		if [[ $all_files == $part_of_list ]];
		then
		    file_tag_exists=true
		fi
	    done
	    if [[ $file_tag_exists == "false" ]];
		then
		echo "###################################################################################################################"
		echo "SNUanalyzer::sktree :: WARNING :: Check input list. "$part_of_list" is not a valid input name for latest snuversion"
		echo "###################################################################################################################"
		echo "Check options using 'sktree -l. If your sample is missing run 'sktree -A' "
		exit 1
	    else
		valid_list_samples=$((valid_list_samples+1))
	    fi
	done
	if [[ $valid_list_samples -eq 0 ]];
	    then
	    echo "SNUanalyzer::sktree :: ERROR :: All input list is invalid. Exiting"
	    exit 1
	fi
	
	if [[ $submit_sampletag != "MC" ]];
	    then
	    
	    if [[ $submit_file_tag  != ""  ]];
		then
		 echo "SNUanalyzer::sktree :: ERROR :: Running with 'sktree -i and -list is not valid'"
		exit 1
	    fi
	fi
	list_isdata=true
	if [[ $submit_sampletag == "MC" ]];
	    then
	    list_isdata=false
	fi
	if [[ $list_isdata == "false" ]];
	    then
	    if [[ $submit_sampletag != ""  ]];
		then
		echo "SNUanalyzer::sktree :: ERROR :: Running with 'sktree -list <list> and -S <samplelist> is not valid unless samplelist is data'"
		exit 1
	    fi
	fi
	
fi

if [[ $submit_sampletag  != ""  ]];
    then
    ARGTEST2=$submit_sampletag
    eval teststreams=(\${$ARGTEST2[@]})
    
    counter=${#teststreams[@]} 
	if [[ $counter -eq 0 ]];
	    then
	    echo "SNUanalyzer::sktree :: ERROR :: Input list is invalid: 'sktree -S samples'"
	    echo "SNUanalyzer::sktree :: HELP :: possible inputs are:"
	    for ilist in  ${list_of_sampletags[@]};
	      do
	      eval samplelist=(\${$ilist[@]})
	      if [[ $ilist == "MC" ]];
		  then
		   echo $ilist"  =  ( All availble MC )"
		  else
		  echo $ilist"  =  (${samplelist[*]})"   
	      fi

	    done
	    exit 1
	fi

	list_isdata=true
	if [[ $submit_sampletag == "MC" ]];
	    then
	    list_isdata=false
	fi
	
	if [[ $list_isdata == "false" ]];
	    then
	    if [[ $submit_file_tag  != ""  ]];
		then
		echo "LSNUnalyzer::sktree :: ERROR :: Running with 'sktree -i <samplename> and -S <samplelist> is not valid unless samplelist is data'"
		 exit 1
	    fi
	
	    if [[ $submit_file_list != ""  ]];
		then
		echo "SNUanalyzer::sktree :: ERROR :: Running with 'sktree -list <list> and -S <samplelist> is not valid unless samplelist is data'"
		exit 1
	    fi
	fi

fi


## RUN PARAMETERS
job_cycle="$submit_analyzer_name"

if [[ $changed_submit_version_tag == "false" ]];
    then
    
    if [[ $submit_analyzer_name == *"SKTreeMaker"* ]];
	then
	submit_version_tag=""
    fi
fi



#### HARDCODE the skinput for sktreemakers
if [[ $submit_analyzer_name == "SKTreeMaker" ]];
    then 
    submit_skinput=false
    job_skim="FLATSNU"
    if [[ $set_sktreemaker_debug == "false" ]];
	then
	job_njobs=1000
    else 
	job_njobs=-311
    fi
    if [[ $submit_version_tag == "" ]];
	then 
	echo "When running SKTreeMaker you need to specify the snuversion: -c "
	echo "Options are: "
	for ic in  ${list_of_snuversions[@]};
	  do
	  echo $ic
	done
	exit 1
    fi
fi
if [[ $submit_analyzer_name == "SKTreeMakerNoCut" ]];
    then
    submit_skinput=false
    job_skim="FLATSNU"
    
    if [[ $set_sktreemaker_debug == "false" ]];
        then
	job_njobs=1000
    else
	job_njobs=-311
    fi
    if [[ $submit_version_tag == "" ]];
        then
        echo "When running SKTreeMaker you need to specify the snuversion: -c "
	echo "Options are: "
	for ic in  ${list_of_snuversions[@]};
          do
          echo $ic
        done

        exit 1
    fi
    
fi

if [[ $submit_analyzer_name == "SKTreeMakerDiLep" ]];
    then
    submit_skinput=true
    job_skim="SKTree_LeptonSkim"
    if [[ $set_sktreemaker_debug == "false" ]];
	then
        job_njobs=1000
    else
	job_njobs=-311
    fi
    if [[ $submit_version_tag == "" ]];
	then
	echo "When running SKTreeMaker you need to specify the snuversion: -c "
	echo "Options are: "
	for ic in  ${list_of_snuversions[@]};
	  do
	  echo $ic
	done
	
	exit 1
    fi
    
fi

if [[ $submit_analyzer_name == "SKTreeMakerTriLep" ]];
    then
    submit_skinput=true
    job_skim="SKTree_LeptonSkim"
    if [[ $set_sktreemaker_debug == "false" ]];
        then
        job_njobs=1000
    else
	job_njobs=-311
    fi
    if [[ $submit_version_tag == "" ]];
	then
	echo "When running SKTreeMaker you need to specify the snuversion: -c "
        echo "Options are: "
        for ic in  ${list_of_snuversions[@]};
          do
          echo $ic
        done

        exit 1
    fi

fi
if [[ $submit_analyzer_name == "SKTreeMakerSS" ]];
    then
    submit_skinput=true
    job_skim="SKTree_LeptonSkim"
    if [[ $set_sktreemaker_debug == "false" ]];
        then
        job_njobs=1000
    else
        job_njobs=-311
    fi
    if [[ $submit_version_tag == "" ]];
	then
        echo "When running SKTreeMaker you need to specify the snuversion: -c "
        echo "Options are: "
        for ic in  ${list_of_snuversions[@]};
          do
          echo $ic
        done

	exit 1
    fi

fi


outdir="/data2/SNU_SKTreeOutput/JobOutPut/"$USER"/SNUanalyzer/"
if [[ ! -d "${outdir}" ]]; then
    mkdir "/data2/SNU_SKTreeOutput/JobOutPut/"$USER
    mkdir ${outdir}
fi


outputdir_snu=$outdir"/data/output/SNU/"

outputdir_analyzer=$outdir"/data/output/SNU/"$submit_analyzer_name
dir_tag="period"$snudatatag"/"

if [[  $job_cycle != "SKTreeMaker"* ]];
    then
    if [[ ! -d "${outputdir_snu}" ]]; then
	mkdir $outdir"/data/"
	mkdir $outdir"/data/output/"
	mkdir ${outputdir_snu}
    fi
    if [[ ! -d "${outputdir_analyzer}" ]]; then
	mkdir ${outputdir_analyzer}
    fi
    if [[ $job_data_lumi  == $snudatatag ]]
	then
	dir_tag="period"$snudatatag"/"
    fi
fi

outputdir_mc=${outputdir_analyzer}"/"${dir_tag}


if [[ $USER == "jalmond" ]];
    then
    
    if [[  $job_cycle != "SKTreeMaker"* ]];
    then
	
	if [[ ! -d "${outputdir_mc}" ]]; then
            mkdir ${outputdir_mc}
            echo "Making " + ${outputdir_mc}
	fi
    fi
    DATE=`date -dlast-monday +%Y-%m-%d`

    if [[ $submit_dir_tag != "NULL" ]]; then
	dir_tag=${dir_tag}$submit_dir_tag"/"
    else
	dir_tag=${dir_tag}${DATE}"/"
    fi
    outputdir_mc=${outputdir_analyzer}"/"${dir_tag}


    if [[  $job_cycle != "SKTreeMaker"* ]];
    then
	
        if [[ ! -d "${outputdir_mc}" ]]; then
            mkdir ${outputdir_mc}
            echo "Making " + ${outputdir_mc}
        fi
    fi
 
    if [[  $job_cycle != "SKTreeMaker"* ]];
    then

        if [[ ! -d "${outputdir_mc}" ]]; then
            mkdir ${outputdir_mc}
            echo "Making " + ${outputdir_mc}
        fi
	echo "Copying "
	cp -r $ANALYZER_DIR/SNUAnalysis/ $outputdir_mc
    fi
fi

outputdir_data=${outputdir_mc}"Data/"
outputdir_np=${outputdir_mc}"Fake/"

if [[ $runDATA == "true" ]];
    then
    if [[ $job_skim == *"No"* ]];
	then echo "SNUanalyzer::sktree :: ERROR :: There are no NoCut skims for data"
	exit 1
    fi
    
    ARGDATASKIM=$submit_sampletag
    eval out_streams_skimcheck=(\${$ARGDATASKIM[@]})
    if [[ $job_skim == *"DiLep"* ]];
    then
	for stream_check in ${out_streams_skimcheck[@]};
	do
	    if [[ $stream_check == *"Single"* ]];
	    then
		echo "SNUanalyzer::sktree :: WARNING :: There are only SPECIAL DiLepton skims for "$stream_check" in snuversion "$submit_version_tag
		#exit 1
	    fi
	done
    fi
    if [[ $job_skim == *"TriLep"* ]];
        then
        for stream_check in ${out_streams_skimcheck[@]};
          do
          if [[ $stream_check == *"Single"* ]];
              then
              echo "SNUanalyzer::sktree :: WARNING :: There are only special TriLepton skims for "$stream_check" in snuversion "$submit_version_tag

              #exit 1
          fi
        done
    fi
    if [[ $job_skim == *"SS"* ]];
    then
        for stream_check in ${out_streams_skimcheck[@]};
          do
          if [[ $stream_check == *"Single"* ]];
              then
              echo "SNUanalyzer::sktree :: WARNING :: There are only special SSepton skims for "$stream_check" in snuversion "$submit_version_tag

              #exit 1                                                                                                                                                                
          fi
	  done
    fi


    logger+=${out_streams_skimcheck[*]}":"
   
fi

workdirtag="/data2/"
if [[ $job_output_dir == "" ]]
    then
    if [[ $runDATA  == "true" ]];
	then
	workdirtag="/data7/DATA/"

	job_output_dir=$outputdir_data
    fi
    if [[ $runMC  == "true" ]];
	then
	job_output_dir=$outputdir_mc
	
    fi
    if [[ $job_run_fake == "True" ]];
	then
	job_output_dir=$outputdir_np
    fi

fi
output_datafile=${outputdir_mc}"/"$job_cycle"_data_snu_"${submit_version_tag}".root"

if [[  $job_cycle != "SKTreeMaker"* ]];
    then

    if [[ ! -d "${outputdir_mc}" ]]; then
	mkdir ${outputdir_mc}
	echo "Making " + ${outputdir_mc}
    fi
    
    if [[ ! -d "${outputdir_data}" ]]; then
	mkdir ${outputdir_data}
	echo "Making " + ${outputdir_data}
    fi
    if [[ ! -d "${outputdir_np}" ]]; then
        mkdir ${outputdir_np}
        echo "Making " + ${outputdir_np}
    fi
    if [[ ! -d "${job_output_dir}" ]]; then
	mkdir ${job_output_dir}
	echo "Making " + ${job_output_dir}
	if [[ ! -d "${job_output_dir}" ]]; then
	        echo "SNUanalyzer::sktree :: ERROR :: Output directory is set by user. "
		    echo "SNUanalyzer::sktree :: ERROR :: Make Directory ried but failed"
    echo "SNUanalyzer::sktree :: ERROR :: User needs to make directory."
fi
    fi
    
fi

if [[ $submit_analyzer_name ==  "" ]];
    then
    echo "No analyzer set: set with -a"
    echo "List of Classes are:"
    while read line
      do
      if [[ $line == *"C++ class"* ]];
          then
          if [[ $line != *"AnalyzerCore"* ]];
              then
              sline=$(echo $line | head -n1 | awk '{print $5}')
              suffix="+;"
              sline=${sline%$suffix}
              echo $sline
          fi
      fi
    done < ${ANALYZER_DIR}/$linkdef_filepath
    exit 1
fi


echo "analyzer::sktree :: INFO :: Analyzer = "${job_cycle}

if [[ $job_run_flip == "True" ]];
    then
    echo "analyzer::sktree :: INFO :: RunFlip = True"
fi

if [[ $job_run_fake == "True" ]];
    then
    echo "analyzer::sktree :: INFO ::RunFake= True"
fi

if [[ $runDATA  == "true" ]];
    then
    ARGPERIOD=data_periods
    eval array_dataperiods=(\${$ARGPERIOD[@]})
    counter=${#array_dataperiods[@]}
    if [[ $counter -eq 0 ]];
            then
	    echo "Number of periods == 0"
	        echo "Problem with initialising 'sktree -p'" 
            exit 1
        fi
fi

if [[ $changed_skim == "true" ]];
    then
    
    if [[ $changed_skinput  == "true" ]];
	then
	
	if [[ $submit_skinput == "false" ]];
	        then
	        if [[ $job_skim != *"SNU"* ]];
		    then 
		    echo "analyzer::sktree :: ERROR :: Invalid skim -s <skim> is set when -sktree <false> is also set. run 'sktree -h' for options "
		    exit 1
		        fi
		else
	        if [[ $job_skim == *"SNU"* ]];
		    then
		    echo "analyzer::sktree :: ERROR :: Invalid skim -s <skim> is set when -sktree <true> is also set. run 'sktree -h' for options "
		    exit 1
		        fi
		    
		fi
	
    fi
    if [[ $changed_skinput  != "true" ]];
	then
	if [[ $job_skim == "FLATSNU" ]];
	        then
	        submit_skinput=false
		else
	        submit_skinput=true
		fi
    fi
    
fi



skim_output_message="analyzer::sktree :: INFO :: Skim set to "${job_skim}" (Default)"

if [[ $changed_skim == "true" ]];
    then
    skim_output_message="analyzer::sktree :: INFO :: Skim set to "${job_skim}
fi

if [[ $job_skim == "FLATSNU" ]];
    then
    echo $skim_output_message
    if [[ $submit_skinput == "true" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set to FLATSNU: Yet you set -sktree <true>"
	echo "Fix submission"
	exit 1
    fi
    
elif [[ $job_skim == "SKTree_NoSkim" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set  to SKTree_NoSkim: Yet you set -sktree <false>"
	echo "Fix submission"
	exit 1
    fi
    
elif [[ $job_skim == "SKTree_LeptonSkim" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set  to SKTree_LeptonSkim: Yet you set -sktree <false>"
echo "Fix submission"
exit 1
    fi
    
elif [[ $job_skim == "SKTree_DiLepSkim" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set  to SKTree_DiLepSkim: Yet you set -sktree <false>"
	echo "Fix submission"
	exit 1
    fi

elif [[ $job_skim == "SKTree_TriLepSkim" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
	then
        echo "analyzer::sktree :: ERROR :: skim set  to SKTree_TriLepSkim: Yet you set -sktree <false>"
	echo "Fix submission"
        exit 1
    fi
elif [[ $job_skim == "SKTree_SSLepSkim" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
        then
        echo "analyzer::sktree :: ERROR :: skim set  to SKTree_SSLepSkim: Yet you set -sktree <false>"
        echo "Fix submission"
	exit 1
    fi


elif [[ $job_skim == "NoCut" ]]
    then
    echo $skim_output_message

    if [[ $submit_skinput == "false" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set  to NoCut: Yet you set -sktree <false>"
	echo "Fix submission"
	exit 1
    fi
    
    
elif [[ $job_skim == "Lepton" ]]
    then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
	then
	echo "analyzer::sktree :: ERROR :: skim set  to Lepton: Yet you set -sktree <false>"
echo "Fix submission"
exit 1
    fi
    
elif [[ $job_skim == "DiLep" ]]
then
    echo $skim_output_message
    if [[ $submit_skinput == "false" ]];
    then
	echo "analyzer::sktree :: ERROR :: skim set  to DiLep: Yet you set -sktree <false>"
	echo "Fix submission"
	exit 1
    fi
    
fi
  
data_lumi_output_message="analyzer::sktree :: INFO :: Running on MC: target lumi set to 13 TeV periods  "${job_data_lumi}" (Default)"
if [[ $chanhed_job_data_lumi == "true" ]];
    then
    data_lumi_output_message="analyzer::sktree :: INFO :: Running on MC: target lumi set to 13 TeV periods  "${job_data_lumi} 
fi

if [[ $runMC  == "true" ]];
    then
    job_data_lumi_exists=false
    for dataperiod in  ${snudataperiods[@]};
    do
    if [[ $job_data_lumi == $dataperiod ]];
        then
	job_data_lumi_exists=true
        if [[ $job_cycle != *"SKTreeMaker"* ]];
        then
            echo $data_lumi_output_message
        fi
    fi
    done
    
    if [[ $job_data_lumi == $snudatatag ]]
	then
	job_data_lumi_exists=true
	if [[ $job_cycle != *"SKTreeMaker"* ]];
	    then
	    echo $data_lumi_output_message
        fi
    elif [[ $job_data_lumi == "GH" ]]
        then
	job_data_lumi_exists=true
	
        if [[ $job_cycle != *"SKTreeMaker"* ]];
            then
            echo $data_lumi_output_message
        fi
    elif [[ $job_data_lumi == "H" ]]
        then
	job_data_lumi_exists=true

        if [[ $job_cycle != *"SKTreeMaker"* ]];
            then
            echo $data_lumi_output_message
        fi

    elif [[ $job_data_lumi == "ALL" ]]
	then
	job_data_lumi_exists=true

	if [[ $job_cycle != *"SKTreeMaker"* ]];
	    then
	    echo $data_lumi_output_message
        fi
    else
	if [[ $job_data_lumi_exists == "false" ]];
	then
	    echo $snudatatag
	    echo "analyzer::sktree :: ERROR :: Invalid setting for data lumi: sktree -p"
	    echo "'sktree -h ' for help."
            exit 1
	fi
    fi
fi

skinput_output_message="analyzer::sktree :: INFO :: SKTree input set to "${submit_skinput}
if [[ $submit_skinput == "true" ]];
    then
    echo $skinput_output_message
    
elif [[ $submit_skinput == "True" ]]
    then
    echo $skinput_output_message
elif [[ $submit_skinput == "false" ]]
        then
    echo $skinput_output_message
elif [[ $submit_skinput == "False" ]]
    then
    echo $skinput_output_message
else
    echo "analyzer::sktree :: ERROR :: Invalid setting for sktree_input: sktree -sktree (run 'sktree -h' for options)"
    exit 1
fi

if [[ $job_loglevel != "ERROR" ]]
    then
    if [[ $job_loglevel != "WARNING" ]]
	then
	if [[ $job_loglevel != "INFO" ]]
	then
	    if [[ $job_loglevel != "DEBUG" ]]
            then
		echo "analyzer::sktree :: ERROR :: loglevel = "${job_loglevel}
		echo "This is invalid: 'sktree -h' for options"
		exit 1
	    fi
	fi
    fi
fi


SNUDEBUG="False"
if [[ $job_loglevel ==  "DEBUG" ]]
then
    SNUDEBUG="True"
fi

#njobs_output_message="analyzer::sktree :: INFO :: Number of subjobs = "${job_njobs}" (Default)"

if [[ $changed_job_njobs == "true" ]];
    then
    
    if [[ $job_cycle != *"SKTreeMaker"* ]];
    then
	if [[ $job_njobs -eq 311 ]];
	    then
	    job_njobs=-300
	elif [[ $job_njobs -gt 200 ]];
	    then
	    if [ $HOSTNAME != "cmscluster.snu.ac.kr" ];
	    then
		echo "analyzer::sktree :: WARNING :: njobs set set out of range (> 200)"
		job_njobs=300
	    fi
	fi
    fi
    njobs_output_message="analyzer::sktree :: INFO :: Number of subjobs = "$job_njobs 
fi

echo $njobs_output_message

loglevel_output_message="analyzer::sktree :: INFO :: loglevel = "${job_loglevel}" (Default)"
if [[ $changed_job_loglevel == "true" ]];
    then
    loglevel_output_message="analyzer::sktree :: INFO :: loglevel = "$job_loglevel
fi
echo $loglevel_output_message


if [[ $submit_version_tag != "" ]];
    then
    job_allowed_snuversion=false
    for job_iac in  ${list_of_snuversions[@]};
      do
      if [[ $job_iac == $submit_version_tag ]];
          then
          job_allowed_snuversion=true
      fi
    done
    if [[ $job_allowed_snuversion == "false" ]];
        then
        echo "analyzer::sktree :: ERROR :: Invalid entry for snuversion"
        echo "'sktree -c <version>'. Options are:"
        for job_iac in  ${list_of_snuversions[@]};
          do
          echo $job_iac
	done
        exit 1
    else 
	if [[ $submit_skinput == "true" ]];
	    then
	    echo "analyzer::sktree :: INFO :: Snuversion for sktrees = "${submit_version_tag}
	elif [[ $submit_skinput == "false" ]]
	    then
	    echo "analyzer::sktree :: INFO :: Snuversion for flat snuuples  = "${submit_version_tag}
	fi
    fi
    else  echo "analyzer::sktree :: INFO :: User opted to check all snuversions for input" 
fi



if [[ $runMC  == "true" ]];
    then

    if [[ $submit_file_tag  != ""  ]];
      then
        echo "analyzer::sktree :: INFO :: Single File to process = "${submit_file_tag}  
	logger+=${submit_file_tag}":"
    fi
    if [[ $submit_file_list  != ""  ]];
	then
        ARGOUT=$submit_file_list
        eval out_input_samples=(\${$ARGOUT[@]})
	echo "analyzer::sktree :: INFO :: ARRAY of samples  to process = "${out_input_samples[*]}
	logger+=${out_input_samples[*]}":"
	    
    fi
fi

if [[ $job_nevents != -1 ]];
    then 
    echo "analyzer::sktree :: INFO :: Set number of events per subjob to : "$job_nevents 
fi

if [[ $job_nskip != -1 ]];
    then
    echo "analyzer::sktree :: INFO :: Set number of events to skip: "$job_nskip
    echo "analyzer::sktree :: INFO :: This is only useful for debugging failed jobs"
fi

if [[ $runDATA  == "true" ]];
    then
    if [[ $submit_sampletag  != ""  ]];
	then
	ARGDATA=$submit_sampletag
	eval out_streams=(\${$ARGDATA[@]})
	echo "analyzer::sktree :: INFO :: ARRAY of samples  to process = "${out_streams[*]}
    fi
    
    ARGPERIOD=data_periods
    eval array_dataperiods=(\${$ARGPERIOD[@]})
    echo "analyzer::sktree :: INFO :: Periods to be processed = "${array_dataperiods[*]}
	
fi

outputdir_output_message="analyzer::sktree :: INFO :: output directory = "${job_output_dir}" (Default)"
 
if [[ $changed_job_output_dir == "true" ]];
    then
    outputdir_output_message="analyzer::sktree :: INFO :: output directory = "${job_output_dir}

fi
if [[ $job_cycle != *"SKTreeMaker"* ]];
    then
    echo $outputdir_output_message
fi


########## GET EMAIL

snu_email=""
while read line
do
    prefix="email = "
    if [[ $line == $prefix* ]];
    then
        line=${line:${#prefix}}
        snu_email=$line
    fi
done < ${ANALYZER_DIR}/bin/snuconfig

if [[ $submit_analyzer_name == *"SKTreeMaker"* ]];
    then
    
    echo "------------------------------------------------------------------------------------------------------" >> sktree_logger.txt
    echo "--User =" $USER " :" >> sktree_logger.txt
    echo "--Date:" >> sktree_logger.txt 
    date >> sktree_logger.txt
    echo "--Comment on reason for making sktrees (first time, bug fix, new variable?):" >> sktree_logger.txt
    echo "-- If testing just add 'TEST' after Message: and this will leave no comment in the change log">> sktree_logger.txt
    echo -e "" >> sktree_logger.txt
    echo "Message:" >> sktree_logger.txt
    if [[ $submit_sk_message != "" ]];	then
	echo $submit_sk_message >> sktree_logger.txt 
    fi
    echo -e "" >> sktree_logger.txt
    echo "--Samples:" >> sktree_logger.txt
    echo $logger >> sktree_logger.txt
    echo "------------------------------------------------------------------------------------------------------" >> sktree_logger.txt
    cp sktree_logger.txt sktree_logger_tmp.txt
    
    snu_editor=""
    while read line
    do

	prefix="editor = "
	if [[ $line == $prefix* ]];
	then
            line=${line:${#prefix}}
	    snu_editor=$line
        fi
    done < ${ANALYZER_DIR}/bin/snuconfig
    echo $snu_editor' sktree_logger.txt'  >> edit.sh
    
    
    if [[ $submit_sk_message == "" ]];
        then
	source edit.sh
	
	if diff sktree_logger.txt sktree_logger_tmp.txt  >/dev/null ; then
	    echo "No comment added: exiting process"
	    rm sktree_logger.txt
	    rm sktree_logger_tmp.txt
	    exit 1
	else
	    echo "Commented added to log:"
	fi
    fi
    
    makelog="True"
    while read line
      do
      if [[ $line == "Message: TEST" ]];
	  then
	  makelog="False"
      fi
      if [[ $line == "Message:TEST" ]];
	  then
	  makelog="False"
      fi
      
    done < sktree_logger.txt
    
    
    if [[ $makelog == "True" ]];
    then
	if [[ -f $ANALYZER_SKTreeLOG_DIR"/$submit_analyzer_name"_${submit_version_tag}".log" ]]; then
	    while read line
	    do
		echo $line >> sktree_logger.txt
	    done < $ANALYZER_SKTreeLOG_DIR"/"$submit_analyzer_name"_${submit_version_tag}.log"
	    
	    echo -e "\n" >> sktree_logger.txt
	    
	    cp sktree_logger.txt $ANALYZER_SKTreeLOG_DIR"/"$submit_analyzer_name"_"${submit_version_tag}".log"
	else
	    cp sktree_logger.txt $ANALYZER_SKTreeLOG_DIR"/"$submit_analyzer_name"_"${submit_version_tag}".log"
	fi
    fi
    rm edit.sh
    rm sktree_logger.txt
    rm sktree_logger_tmp.txt
fi

if [[ $GetOutPutDir == "True" ]];
then
    if [[ $submit_analyzer_name ==  "" ]];
    then
	echo "No analyzer set: set with -a"
	echo "List of Classes are:"
    while read line
    do
	if [[ $line == *"C++ class"* ]];
        then
            if [[ $line != *"AnalyzerCore"* ]];
            then
		sline=$(echo $line | head -n1 | awk '{print $5}')
		suffix="+;"
		sline=${sline%$suffix}
		echo $sline
            fi
	fi
    done < ${ANALYZER_DIR}/$linkdef_filepath
    exit 1
    fi
    
    echo "Outputdir for "${submit_analyzer_name}" is "$outputdir
    exit 1
fi

#iddir=/data2/SNU_SKTreeOutput/$USER/JobID/ >> declared in snuconf now
if [[ ! -d "${iddir}" ]]; then
    mkdir $iddir
fi

qnameisok="false"
for qname in ${queueoptions[@]};
do
    if [[ $queuename  == *${qname}* ]]; then 
	qnameisok="true"
    fi

done
if [[ $qnameisok == "false" ]];then
    echo "--------------------------------------------------------------------------"
   
    echo "Error in Queuename: run  sktree -qlist"
    echo "--------------------------------------------------------------------------"
    exit 1
fi

runboth="false"
if [[ $runDATA  == "true" ]];
    then
    if [[ $runMC  == "true" ]];
    then
	runboth="true"
	echo "analyzer::sktree :: INFO :: runing MC+DATA, DATA will be run in background"
    fi
fi


if [[ $runDATA  == "true" ]];
    then
    

    ARG1=$submit_sampletag
    eval streams=(\${$ARG1[@]})
    for istream in ${streams[@]};
      do
      echo "analyzer::sktree :: INFO :: stream = "${istream}

      if [[ $istream == "" ]];
	  then
	  continue
      fi


      JobID=""
      if [ -f $JOBIDTAG ];
      then
	  cp $JOBIDTAG $iddir/IDtmp.txt
	  
	  while read line
	  do
              sline1=$(echo $line | head -n1 | awk '{print $1}')
              sline2=$(echo $line | head -n1 | awk '{print $2}')
              sline2=$((sline2+1))
              new_line=$sline1" "$sline2
              echo $new_line > $JOBIDTAG
              JobID=$sline2
	  done < $iddir/IDtmp.txt
	  
      else
	  if [[ $USER == "jalmond" ]];
	  then
              echo "Error in  $JOBIDTAG"
              exit 1
	  fi
	  JobID=$RANDOM
      fi

      if [ "$JobID" -eq "$JobID" ] 2>/dev/null; then
          tagger=$JobID
      else
          tagger=$RANDOM
      fi
      
      
      ### set all inputs to default in functions.sh
      source functions.sh


      run_in_bkg=${job_run_bkg}
      #### change input to those specified in sktree command
      if [[ $runboth  == "false" ]];
      then
	  run_in_bkg=${job_run_bkg}
	  
      else
	  sendemail="True"
      fi
      echo "analyzer::sktree :: INFO :: queue ="$queuename
      queue=${queuename}
      nevents=${job_nevents}
      skipevent=${job_nskip}
      runnp=${job_run_fake}
      runtau=${job_run_taus}
      runcf=${job_run_flip}
      cycle=${job_cycle}
      skinput=${submit_skinput}
      useskim=${job_skim}
      skflag=${submit_skflag}
      njobs=$job_njobs
      data_lumi=$job_data_lumi
      loglevel=$job_loglevel
      logstep=$job_logstep
      stream=${istream}
      drawhists=${submit_draw}
      setnjobs=${changed_job_njobs}
      submitallfiles=${job_submitallfiles}
      tmpfilename=${job_tmp_filename}
      if [[ $changed_job_output_dir == "true" ]];
	  then
	  outputdir=${job_output_dir}
      else
	  if [[ $job_run_fake == "True" ]];
	      then
	      outputdir=$outputdir_np
	      else 
	      outputdir=$outputdir_data
	  fi
      fi
      snuversion=${submit_version_tag}
    
      ARG=data_periods
      eval input_samples=(\${$ARG[@]})

      
      
      echo "analyzer::sktree :: INFO :: JOB ID["$tagger"]"
      if [[ ${run_in_bkg} == "True" ]];
	  then
	  echo "Job sent to $workdirtag/SNU_SKTreeOutput/"$USER"/CLUSTERLOG"$tagger"/"${tagger}
	  echo "Job summary will be sent via email when complete:"
	  logdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/
	  if [[ ! -d "${logdir}" ]]; then
	      mkdir $logdir
	  fi
	  joblogdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/${tagger}/
	  if [[ ! -d "${joblogdir}" ]]; then
	      mkdir $joblogdir
	  fi
	  source submit.sh $tagger &>  $joblogdir/job.txt &

      else
	  source submit.sh $tagger
      fi
      
      if [[ ${job_run_fake} == "True" ]];
      then
	  if [[ ${run_in_bkg} != "True" ]];then
	      mergefake $istream
	  fi
      else
	  if [[ ${run_in_bkg} != "True" ]];then
              mergeoutput $istream
	  fi
      fi
    done
fi

if [[ $runMC  == "true" ]];

    then

    JobID=""
    
    if [ -f $JOBIDTAG ];
    then
	cp $JOBIDTAG $iddir/IDtmp.txt
	
	while read line
	do
            sline1=$(echo $line | head -n1 | awk '{print $1}')
            sline2=$(echo $line | head -n1 | awk '{print $2}')
	    sline2=$((sline2+1))
            new_line=$sline1" "$sline2
            echo $new_line > $JOBIDTAG
            JobID=$sline2
	done < $iddir/IDtmp.txt
	
    else
	if [[ $USER == "jalmond" ]];
	then
            echo "Error in  $JOBIDTAG"
            exit 1
	fi
	JobID=$RANDOM
    fi

    
    tagger=""
    if [ "$JobID" -eq "$JobID" ] 2>/dev/null; then
        tagger=$JobID
    else
        tagger=$RANDOM
    fi
    
    
    ### set all inputs to default in functions.sh
    source functions.sh
    
    if [[ $runboth  == "false" ]];
      then
        run_in_bkg=${job_run_bkg}
      else
        sendemail="True"
    fi


    #### change input to those specified in sktree command
    run_in_bkg=${job_run_bkg}
    cycle=${job_cycle}
    echo "analyzer::sktree :: INFO :: queue ="$queuename
    queue=${queuename}
    skinput=${submit_skinput}
    useskim=${job_skim}
    njobs=$job_njobs
    setnjobs=${changed_job_njobs}
    skflag=${submit_skflag}
    data_lumi=$job_data_lumi
    loglevel=$job_loglevel
    logstep=$job_logstep
    outputdir=${job_output_dir}
    snuversion=${submit_version_tag}
    nevents=${job_nevents}
    skipevent=${job_nskip}
    runnp=${job_run_fake}
    runtau=${job_run_taus}
    runcf=${job_run_flip}
    drawhists=${submit_draw}
    submitallfiles=${job_submitallfiles}
    tmpfilename=${job_tmp_filename}



    if [[ $submit_file_tag  != ""  ]];
        then
	
        declare -a input_samples=("${submit_file_tag}")
	echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

	echo "analyzer::sktree :: INFO :: JOB MC ID["$tagger"]"


	if [[ ${run_in_bkg} == "True" ]];
        then
	    echo "Job sent to $workdirtag/SNU_SKTreeOutput/"$USER"/CLUSTERLOG"$tagger"/"${tagger}
            echo "Job summary will be sent via email when complete:"
	    logdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/
            if [[ ! -d "${logdir}" ]]; then
		mkdir $logdir
	    fi
            joblogdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/${tagger}/
            if [[ ! -d "${joblogdir}" ]]; then
		mkdir $joblogdir
            fi
            source submit.sh $tagger &>  $joblogdir/job.txt &

	else
            source submit.sh $tagger
	fi
    fi
    
    if [[ $submit_file_list  != ""  ]];
	then
	ARG=$submit_file_list
	eval input_samples=(\${$ARG[@]})

	echo "analyzer::sktree :: INFO :: MC JOB ID["$tagger"]"

	if [[ ${run_in_bkg} == "True" ]];
        then
            echo "Job sent to $workdirtag/SNU_SKTreeOutput/"$USER"/CLUSTERLOG"$tagger"/"${tagger}
            echo "Job summary will be sent via email when complete:"
	    logdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/
            if [[ ! -d "${logdir}" ]]; then
		mkdir $logdir
	    fi
            joblogdir=$workdirtag/SNU_SKTreeOutput/$USER/CLUSTERLOG$tagger/${tagger}/
            if [[ ! -d "${joblogdir}" ]]; then
		mkdir $joblogdir
            fi
            source submit.sh $tagger &>  $joblogdir/job.txt &
	    	    
	else
            source submit.sh $tagger
	fi
    fi
fi

