#!/bin/sh
### sets all configurable variables to defaul values

function usage
{

    echo "usage: sktree [-a analyzer] [-S samples] [-i input_file ]"
    echo "              [-s skim] [-list file_array] [-p data_period] [-q queue]"
    echo "              [-d debug_mode] [-c snuversion] [-o outputdir] [-qlist] "
    echo "              [-events number of events] [-nskip events_to_skip] [-ac allversion] [-b run_in_bkg]"
    echo "              [-fake runfake ] [-flip runflip] [-attachhist drawhist]"     
    echo "              [-h (more/debug)][-l <args> ][-g <args>] [-A <args>]"
    echo "              [-D <snuversion>] [-miniaod input_file ] [-xsec input_file] [-efflumi input_file] [-userflag flag]"
    echo "              [-tagdiff <tagname>  -sktreelog  -printID IDNAME -updateselection <object>]   "
    echo "              [-filename xxx]"
}



function rungroupedlist
{
    echo ""
    runFull=false
    if [[ $submit_searchlist == "list" ]];
	then
	runFull=true
    elif [[ $submit_searchlist != "" ]];
	then
	echo "Invalid option for 'sktree -g': Use 'sktree -g' or 'sktree -g list' "
	exit 1
    fi
    
    while read line
      do
      if [[ $line == *"declare"* ]];
	  then
	  sline=$(echo $line | cut -d "" -f2-)
	  if [[ $runFull == "false" ]];
	      then
	      sline=${sline/"="/" "}
	      sline2=$(echo $sline | head -n1 | awk '{print $3}')
	      echo $sline2
	  else 
	      echo $sline
	      echo ""
	  fi
	  
      fi
    done < ${ANALYZER_DIR}/SNURun/txt/list_all_mc_${SNUVERSION}.sh
    
    echo ""
    echo "Arrays made by user"

    while read line
      do
      if [[ $line == *"declare"* ]];
          then
          sline=$(echo $line | cut -d "" -f2-)
	  if [[ $runFull == "false" ]];
	      then
	      sline=${sline/"="/" "}
	      sline2=$(echo $sline | head -n1 | awk '{print $3}')
	      echo $sline2
	      
	  else
	      echo $sline
	      echo ""
          fi
      fi
    done < ${ANALYZER_DIR}/SNURun/txt/list_user_mc.sh

}

getinfo_tag=0
getinfo_string=""


#function checkdata
#{
#    if [[ $submit_sampletag == "" ]];
#	then
#	echo "ERROR in checkdata"
#	exit 1
#    fi
#    eval check_list=(\${$submit_sampletag[@]})
#	
#    for idlist in  ${check_list[@]};##
#
#      do
#      data_isok=false
#      slined="/data2/DATA/cattoflat/Data/"
#      if [[ $job_skim == *"No"* ]];
#          then
#          echo "SNUanalyzer::sktree :: ERROR :: There are no NoCut skims for "$idlist
#          echo  "SNUanalyzer::sktree :: HELP :: Change skim"
#          exit 1
#      fi
#      if [[ $job_skim == *"Lepton"* ]];
#          then
#          slined="/data2/CatNtuples/"${submit_version_tag}"/SKTrees/Data/"
#      fi
#      if [[ $job_skim == *"DiLep"* ]];
#          then
#          slined="/data2/CatNtuples/"${submit_version_tag}"/SKTrees/DataDiLep/"
#      fi
#      
#      dirtag_date=$idlist
#      if [[ $submit_version_tag == *"v7-4" ]];
#	  then
#	  if [[ $idlist == "emu" ]]
#	      then
#	      dirtag_date=ElectronMuon
#	  fi
#	  if [[$idlist== "electron" ]]
#	      then
#	      dirtag_date=DoubleElectron
#	  fi
#	  if [[$idlist== "muon" ]];
#	      then
#	      dirtag_date=DoubleMuon
#	  fi
#	  if [[$idlist== "singlemuon" ]];
#	      then
#	      dirtag_date=SingleMuon
#	  fi
#	  if [[$idlist== "singleelectron" ]];
#	      then
#	      dirtag_date=SingleElectron
#	  fi
#      fi
#      slined=$slined$dirtag_date
#      
#      if [[ -d "${slined}" ]]; then
#          if test "$(ls -A "$slined")"; then
#              data_isok=true
#          fi
#      fi
#      if [[ $data_isok == "false" ]];
#          then
#          echo "Data sample "$dirtag_date" is not available for catversion "$submit_version_tag" and skim "$job_skim
#          exit 1
#      fi
#    done
#}

function getinfo_dataset
{ 
    
    if [[ $submit_file_tag == *$search_tag* ]];
        then
        tmp_submit_file_tag=$submit_snuvlist
        tmp_submit_snuvlist=$submit_file_tag
        submit_snuvlist=$tmp_submit_snuvlist
        submit_file_tag=$tmp_submit_file_tag
    fi
    if [[ $submit_snuvlist != *$search_tag* ]];
        then
        if [[ $submit_snuvlist != "" ]];
            then
            tmp_submit_file_tag=$submit_snuvlist
            tmp_submit_snuvlist=$submit_file_tag
            submit_snuvlist=$tmp_submit_snuvlist
            submit_file_tag=$tmp_submit_file_tag
        fi
    fi

    allowed_snuversion=false
    for iac in  ${list_of_snuversions[@]};
      do
      if [[ $iac == $submit_snuvlist ]];
          then
          allowed_snuversion=true
      fi
    done
    if [[ $submit_snuvlist != "" ]];
        then
        if [[ $allowed_snuversion == "false" ]];
            then
            echo "SNUanalyzer::sktree :: ERROR :: Snuversion "$submit_snuvlist" is not allowed"
	    exit 1
	elif [[ $submit_snuvlist == *"v7-4"*  ]];
            then
            echo "SNUanalyzer::sktree :: ERROR :: 'sktree -D' only works for v7-6-2 and newer"
            exit 1

        fi
    fi

    if [[ $submit_snuvlist  == "" ]];
        then
        submit_snuvlist=${SNUVERSION}
    fi
    if [[ $submit_file_tag == "" ]];
	then
	echo "Need a dataset is input"
	exit 1
    fi
    

    while read line
      do
      if [[ $line == *$getinfo_string* ]];
          then

	  if [[ $line == *" "$submit_file_tag" "* ]];
	      then
	      sline=""
	      if [[ $getinfo_tag == "2" ]];
		  then
		  sline=$(echo $line | head -n1 | awk '{print $2}')
	      elif [[ $getinfo_tag == "3" ]];
                  then
		  sline=$(echo $line | head -n1 | awk '{print $3}')
	      elif [[ $getinfo_tag == "4" ]];
		  then
		  sline=$(echo $line | head -n1 | awk '{print $4}')
	      elif [[ $getinfo_tag == "5" ]];
		  then
                  sline=$(echo $line | head -n1 | awk '{print $5}')
	      else 
		  sline=$(echo $line | head -n1 | awk '{print $6}')
		  
	      fi
	      echo ${sline}
      
	  fi
      fi
      
    done < ${TXTPATH}"SNU_mc_"${submit_snuvlist}".txt"

}



function getalldatasetinfo
{
    allowed_snuversion=false
    for iac in  ${list_of_snuversions[@]};
      do
      if [[ $iac == $submit_snuvlist ]];
          then
          allowed_snuversion=true
      fi
    done
    if [[ $submit_snuvlist != "" ]];
        then
        if [[ $allowed_snuversion == "false" ]];
            then
            echo "SNUanalyzer::sktree :: ERROR :: Snuversion "$submit_snuvlist" is not allowed"
            exit 1
	elif [[ $submit_snuvlist == *"v7-4"*  ]];
            then
            echo "SNUanalyzer::sktree :: ERROR :: 'sktree -D' only works for v7-6-2 and newer"
            exit
        fi
    fi
    
    if [[ $submit_snuvlist  == "" ]];
        then
        submit_snuvlist=${SNUVERSION}
    fi
    echo "######################################################################################"
    echo "SNUAnalyzer::sktree :: INFO :: Looking up information for "$submit_snuvlist 
    echo "######################################################################################"
    echo ""
    echo "https://github.com/vallot/SNUTools/blob/"$submit_snuvlist"/SnuProducer/python/snuDefinitions_cfi.py"
    echo ""

    curl https://github.com/vallot/SNUTools/blob/"$submit_snuvlist"/SnuProducer/python/snuDefinitions_cfi.py >> gitfile.txt
    echo "######################################################################################"


    echo ""
    declare -a ProdInfo=( "globalTag_mc" "globalTag_rd" "bunchCrossing" "lumiJSON" "lumiJSONSilver" "pileupMCmap" "JetEnergyCorrection" )
    
    while read line
      do
      
      for ipr in  ${ProdInfo[@]};
      	do
	if [[  $line == *$ipr* ]];
	  then
	    
	    if [[ $ipr != "bunchCrossing" ]];
		then
		sline2=$(echo $line | head -n1 | awk '{print $5}')
		sline2=${sline2:14}
		
		sline3=$(echo $line | head -n1 | awk '{print $10}')
		sline4=${sline3:27}
		sline5=${sline4//<span/" "}
		if [[ $sline2 == $ipr ]];
		    then
		    echo $ipr"=" $sline5
		fi
	    else
		sline2=$(echo $line | head -n1 | awk '{print $5}')
		sline2=${sline2:14}

                sline3=$(echo $line | head -n1 | awk '{print $9}')
                sline4=${sline3:14}
                sline5=${sline4//td>/" "}
		sline5=${sline5////" "}
                sline5=${sline5//</" "}
		sline5=${sline5//span>/" "}

                if [[ $sline2 == $ipr ]];
                    then
                    echo $ipr"=" $sline5
                fi

		
	    fi
	    
	fi
      done
    done <  gitfile.txt

    rm gitfile.txt

}
function getdatasetname
{
    getinfo_tag=3
    getinfo_string="DATASET:"
    getinfo_dataset
}

function getdatasetxsec
{
    getinfo_tag=4
    getinfo_string=$FLATSNU_MC
    getinfo_dataset
}
function getdatasetefflumi
{
    getinfo_tag=5
    getinfo_string=$FLATSNU_MC

    getinfo_dataset
}

function print_tag_diff
{
    if [[ ${submit_snu_tag} == *$search_tag* ]];then
	
	if [[ ${submit_snu_tag2} == *$search_tag* ]];then
	    print_tag_diff_twotags 
	else
	    print_tag_diff_vs_currenttag	    
	fi
    fi

}

function print_tag_diff7
{
    if [[ ${submit_snu_tag} == *"v7"* ]];then

        if [[ ${submit_snu_tag2} == *"v7"* ]];then
            print_tag_diff_twotags7
        else
            print_tag_diff_vs_currenttag
        fi
    fi

}



function print_tag_diff_twotags
{
    declare -a NEWTAGS=()
    foundtag=False
    while read line
      do
      if [[ $line == *"$submit_snu_tag2"* ]];
	  then
	  sline=$(echo $line | head -n1 | awk '{print $1}')
	  NEWTAGS+=(${sline})
	  foundtag=True
      else
	  if [[ $foundtag == "False" ]]; then
	      continue
	  else
	      if [[ $line == $submit_snu_tag ]]; then
		  break
	      fi
	      NEWTAGS+=(${line})
	  fi
      fi
    done < ${SNUTAGDIR}/LatestTag94X.txt

    for ntag in  ${NEWTAGS[@]};
      do
      echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
      echo "Tag: " $ntag  "(summary of changes wrt previous tag)"
      echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
      
      while read line
	do
	echo $line
      done <  ${SNUTAGDIR}/TagDiff_${ntag}.txt
    done
    
}


function print_tag_diff_twotags7
{
    declare -a NEWTAGS=()
    foundtag=False
    while read line
      do
      if [[ $line == *"$submit_snu_tag2"* ]];
          then
          sline=$(echo $line | head -n1 | awk '{print $1}')
          NEWTAGS+=(${sline})
          foundtag=True
      else
          if [[ $foundtag == "False" ]]; then
              continue
          else
              if [[ $line == $submit_snu_tag ]]; then
                  break
              fi
              NEWTAGS+=(${line})
          fi
      fi
    done < ${SNUTAGDIR}/LatestTag.txt

    for ntag in  ${NEWTAGS[@]};
      do
      echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
      echo "Tag: " $ntag  "(summary of changes wrt previous tag)"
      echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

      while read line
        do
        echo $line
      done <  ${SNUTAGDIR}/TagDiff_${ntag}.txt
    done

}


function update_selection
{
    python ${ANALYZER_DIR}/python/UpdateSelection.py -s $object
}
function printid
{
    
    if [[ $idname == "" ]];then
	echo "No input given:"
    else
	if [[ $idname2 == "" ]];
	then
	    python ${ANALYZER_DIR}/python/PrintIDSelection.py --id1 $idname 
	else
	    python ${ANALYZER_DIR}/python/PrintIDSelection.py --id1 $idname --id2 $idname2
	fi
    fi
    
}


function print_sktreemaker_logfile
{
    if [[ -f /data1/SNUAnalyzer_rootfiles_for_analysis/SNUSKTreeMaker/${submit_analyzer_name}_${submit_snuvlist}.log ]]; then
	while read line
	  do
	  echo $line
	done < /data1/SNUAnalyzer_rootfiles_for_analysis/SNUSKTreeMaker/${submit_analyzer_name}_${submit_snuvlist}.log 
    else
	echo "Invalid input:"
	echo "sktree -sktreelog <analyzername> <snuversion> (i.ie, sktree -sktreelog SKTreeMakerDiLep v7-6-4)"
    fi

}
function print_tag_diff_vs_currenttag
{
    tag_diff_file=$SNUTAGDIR/TagDiff_${submit_snu_tag}.txt
    
    latest_tag=""
    while read line
      do
      if [[ $line == *"HEAD"* ]];
	  then
	  sline=$(echo $line | head -n1 | awk '{print $1}')
	  latest_tag=$sline
      fi
    done < $SNUTAGDIR/LatestTag80X.txt
    
    if [[ $latest_tag == $SNUTAG ]];then
	
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	echo "Input tag is latest tag"
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    else
	
	
	echo "Summary of difference between current tag and "$submit_snu_tag
	declare -a NEWTAGS=()
	while read line
	  do
      if [[ $line == *"HEAD"* ]];
	  then
	  sline=$(echo $line | head -n1 | awk '{print $1}')
	  NEWTAGS+=(${sline})
      else
	  if [[ $line == $SNUTAG ]]; then
	      break
	  fi
	  NEWTAGS+=(${line})
      fi
	done < $SNUTAGDIR/LatestTag801.txt
	
	for ntag in  ${NEWTAGS[@]};
	  do
	  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	  echo "Tag: " $ntag  "(summary of changes wrt previous tag)"
	  echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	  
	  while read line
	    do
	    echo $line
	  done < $SNUTAGDIR/TagDiff_${ntag}.txt
	done
    fi
}

function listavailable
{
    
    if [[ $submit_searchlist == *$search_tag* ]];
        then
        tmp_submit_searchlist=$submit_snuvlist
        tmp_submit_snuvlist=$submit_searchlist
        submit_snuvlist=$tmp_submit_snuvlist
        submit_searchlist=$tmp_submit_searchlist
    fi
    if [[ $submit_snuvlist != *$search_tag* ]];
	then
	if [[ $submit_snuvlist != "" ]];
	    then
	    tmp_submit_searchlist=$submit_snuvlist
	    tmp_submit_snuvlist=$submit_searchlist
	    submit_snuvlist=$tmp_submit_snuvlist
	    submit_searchlist=$tmp_submit_searchlist
	fi
    fi
    
    allowed_snuversion=false
    for iac in  ${list_of_snuversions[@]};
      do
      if [[ $iac == $submit_snuvlist ]];
          then
          allowed_snuversion=true
      fi
    done
    if [[ $submit_snuvlist != "" ]];
	then
        if [[ $allowed_snuversion == "false" ]];
            then
            echo "SNUanalyzer::sktree :: ERROR :: Snuversion "$submit_snuvlist" is not allowed"
            exit 1
        fi
    fi
    
    specified_snuversion=true
    if [[ $submit_snuvlist  == "" ]];
        then
        submit_snuvlist=${SNUVERSION}
        specified_snuversion=false
    fi
    
    echo ""
    echo "List of available samples at SNU. With snuversion " ${submit_snuvlist}
    echo ""
    echo "Samplename  --> datasetname"  
    while read line
      do
      if [[ $line == *$FLATSNU_MC* ]];
	  then
	  if [[ $submit_searchlist == "" ]];
	      then
	      sline1=$(echo $line | head -n1 | awk '{print $1}')
	      sline2=$(echo $line | head -n1 | awk '{print $6}')
	      echo ${sline1} "--> "  ${sline2}
	  fi
	   if [[ $submit_searchlist != "" ]];
              then
	       if [[ $line == *${submit_searchlist}* ]];
		   then
		   sline1=$(echo $line | head -n1 | awk '{print $1}')
		   sline2=$(echo $line | head -n1 | awk '{print $6}')
		   echo ${sline1} "--> "  ${sline2}
	       fi
	   fi
      fi
    done < ${TXTPATH}"SNU_mc_"${submit_snuvlist}".txt"

    echo ""
    echo ""
    echo ""
    echo "List of samples not available in latest available snuversion are:"
    echo ""
    echo "Missing: since miniAOD not available"

    while read line
      do
      if [[ $line == *"Missing:"* ]];
	  then
	  sline=$(echo $line | head -n1 | awk '{print $2}')
	  sline2=$(echo $line | head -n1 | awk '{print $3}')
	  if [[ $submit_searchlist == "" ]];
              then
	      echo ${sline} "--> "  ${sline2}
	  fi
	  if [[ $submit_searchlist != "" ]];
	      then
	      if [[ $line == *${submit_searchlist}* ]];
		  then
		  echo ${sline} "--> "  ${sline2}
	      fi
	  fi
      fi
    done < ${TXTPATH}"SNU_mc_"${submit_snuvlist}".txt"

    echo ""
    echo ""
    echo ""
    echo "Snuuples available at kisti: Can create flatsnuuples->sktrees on request: type 'sktree -r datasetname' "

    while read line
      do
    if [[ $line == *"Available[not produced]"* ]];
	then
	sline=$(echo $line | head -n1 | awk '{print $3}')
	if [[ $submit_searchlist == "" ]];
	    then
	    echo ${sline}
	fi
	if [[ $submit_searchlist != "" ]];
	    then
	    if [[ $line == *${submit_searchlist}* ]];
		then
		echo ${sline}
	    fi
	fi
    fi
    done < ${TXTPATH}"SNU_mc_"${submit_snuvlist}".txt"

    echo ""
    echo ""
    if [[ $specified_snuversion == "true" ]];
	then
	echo "To check availability of SKTrees"
	echo "For No skim run 'sktree -L SKTree_NoSkim " ${submit_searchlist} ${submit_snuvlist} "'" 
	echo "For Lepton skim run 'sktree -L SKTree_LeptonSkim " ${submit_searchlist} ${submit_snuvlist} "'" 
	echo "For DiLepton skim run 'sktree -L SKTree_DiLepSkim " ${submit_searchlist} ${submit_snuvlist} "'" 
	echo "For TriLepton skim run 'sktree -L SKTree_TriLepSkim " ${submit_searchlist} ${submit_snuvlist} "'" 
	echo "For SSLepton skim run 'sktree -L SKTree_SSLepSkim " ${submit_searchlist} ${submit_snuvlist} "'" 
	
    fi
    if [[ $specified_snuversion != "true" ]];
        then
	echo "To check availability of SKTrees"
	echo "For No skim run 'sktree -L SKTree_NoSkim " ${submit_searchlist} "'"
	echo "For Lepton skim run 'sktree -L SKTree_LeptonSkim " ${submit_searchlist}  "'"
	echo "For DiLepton skim run 'sktree -L SKTree_DiLepSkim " ${submit_searchlist}  "'"
	echo "For TriLepton skim run 'sktree -L SKTree_TriLepSkim " ${submit_searchlist}  "'"
	echo "For SSLepton skim run 'sktree -L SKTree_SSLepSkim " ${submit_searchlist}  "'"
	
    fi

    
}



function sendrequestsnu
{
    source mail_snu.sh $request_sample
    snu email.txt | mail -s "SNUTuple request" jalmond@cern.ch
    rm email.txt
}

function listqueue
{
    
    qstat -f 
    echo " " 
    echo "To select a certain queue use -q <qname> like "
    while read line
      do
	if [[ $line == *"###"* ]];then
	    echo ""
	else
	    echo $line
	fi
    done < /data1/SNUAnalyzer_rootfiles_for_analysis/SnutupleConfig/QUEUE/queuelist.txt


}


function sendrequest
{
    source mail.sh $request_sample
    snu email.txt | mail -s "SKTree request" jalmond@cern.ch
    rm email.txt
}

function runlist
{
    
    if [[ $submit_searchlist == *$search_tag* ]];
	then
	tmp_submit_searchlist=$submit_snuvlist
	tmp_submit_snuvlist=$submit_searchlist
	submit_snuvlist=$tmp_submit_snuvlist
	submit_searchlist=$tmp_submit_searchlist
    fi
    if [[ $submit_snuvlist != *$search_tag* ]];
	then
        if [[ $submit_snuvlist != "" ]];
            then
            tmp_submit_searchlist=$submit_snuvlist
            tmp_submit_snuvlist=$submit_searchlist
            submit_snuvlist=$tmp_submit_snuvlist
            submit_searchlist=$tmp_submit_searchlist
        fi
    fi
    
    allowed_snuversion=false
    for iac in  ${list_of_snuversions[@]};
      do
      if [[ $iac == $submit_snuvlist ]];
	  then
	  allowed_snuversion=true
      fi
    done
    if [[ $submit_snuvlist != "" ]];
	then
	if [[ $allowed_snuversion == "false" ]];
	    then
            echo "SNUanalyzer::sktree :: ERROR :: Snuversion "$submit_snuvlist" is not allowed"
	    exit 1
	fi
    fi

    specified_snuversion=true
    if [[ $submit_snuvlist  == "" ]];
	then
	submit_snuvlist=${SNUVERSION}
	specified_snuversion=false
    fi
    
    check_path=""
    
    if [[ $submit_skim  == "" ]];
        then
	echo "sktree -L <skim>"
	echo "Need to set skim: Options are FLATSNU/SKTree_NoSkim/SKTree_LeptonSkim/SKTree_DiLepSkim/SKTree_TriLepSkim"
	echo "Can also specify snuversion AND/OR search filter list (after skim)"
	echo "example 1) sktree -L FLATSNU"
	echo "example 2) sktree -L FLATSNU v7-6-3"
	echo "example 3) sktree -L SKTree_LeptonSkim QCD v7-6-2"
	echo "example 4) sktree -L SKTree_DiLepSkim v7-6-3 DY"
	exit 1
    fi

    if [[ $submit_skim  == "FLATSNU" ]];
    then
        check_path=$FLATSNU_MC
    fi
    
    isNoCut=false
    isLepton=false
    isDiLep=false
    isHNDiLep=false
    isHNFake=false
    isHNFatJet=false
    isTriLep=false
    isSSLep=false
    if [[ $submit_skim  == "SKTree_NoSkim" ]]; then 
	isNoCut=true 
    fi
    if [[ $submit_skim  == "NoSkim" ]];
        then
	isNoCut=true
    fi
    if [[ $submit_skim  == "SKTree_LeptonSkim" ]];
	then
	isLepton=true
    fi
    if [[ $submit_skim  == "Lepton" ]];
	then
	isLepton=true
    fi
    if [[ $submit_skim  == "SKTree_DiLepSkim" ]];
        then
	isDiLep=true
    fi
    if [[ $submit_skim  == "SKTree_TriLepSkim" ]];
        then
        isTriLep=true
    fi
    if [[ $submit_skim  == "SKTree_SSLepSkim" ]];
	then
	isSSLep=true
    fi
    if [[ $submit_skim  == "DiLep" ]];
	then
	isDiLep=true
    fi
    if [[ $isNoCut  == "true" ]];
	then
        check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCNoCut/"
    fi
    if [[ $isLepton  == "true" ]];
	then
	check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MC/"
	if [[ ${submit_snuvlist} == *"v7-4-4"* ]];
            then
            check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/Sep15/MC/"
        fi
    fi
    if [[ $isDiLep  == "true" ]];
	then
	check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCDiLep"
	if [[ ${submit_snuvlist} == *"v7-4-4"* ]];
	    then
	    check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/Sep15/MCDiLep"
	fi
    fi
    if [[ $isHNDiLep  == "true" ]];
     then
        check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCHNDiLep"
        if [[ ${submit_snuvlist} == *"v7-4-4"* ]];
            then
            check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/Sep15/MCHNDiLep"
        fi
    fi
    if [[ $isHNFake  == "true" ]];
     then
        check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCHNFake"

    fi

    if [[ $isHNFatJet  == "true" ]];
     then
        check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCHNFatJet"

    fi

    if [[ $isTriLep  == "true" ]];
        then
        check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCTriLep"
        if [[ ${submit_snuvlist} == *"v7-4-4"* ]];
            then
            check_path=""
        fi
    fi
    if [[ $isSSLep  == "true" ]];
	then
	check_path=$SKTREE_MC${submit_snuvlist}"/SKTrees/MCSS"
        if [[ ${submit_snuvlist} == *"v7-4-4"* ]];
            then
            check_path=""
        fi
    fi
    echo $check_path
    if [[ $check_path == "" ]];
	then
	echo "Invalid option for ntuple version: "
	echo "sktree -L <skim>"
        echo "Need to set ntuple version: Options are FLATSNU/SKTree_NoSkim/SKTree_LeptonSkim/SKTree_DiLepSkim/SKTree_SSLepSkim"
	exit 1
    fi
    

    echo "List of samplenames for skim " $submit_skim " available in snutuple version " $submit_snuvlist " are:"
    
    declare -a LISTOFSAMPLES=()
    declare -a UNPROCESSED=()
    
    while read line
      do
      if [[ $line == *${check_path}* ]];
	  then
	  sline=$(echo $line | head -n1 | awk '{print $1}')
	  sline2=$(echo $line | head -n1 | awk '{print $6}')

	  if [[ $submit_searchlist == "" ]];
	      then
	      prefix="SK"
	      suffix="_dilep"
	      suffixhn="_hndilep"
	      suffixhnfake="_hnfake"
	      suffixhnfatjet="_hnfatjet"
              suffix2="_trilep"
              suffixss="_sslep"

	      if [[ $sline == *${prefix}* ]];
		  then
		  sline=${sline:2}		  
	      fi
	      if [[ $sline == *${suffix}* ]];
                  then
		  sline=${sline%$suffix}
	      fi
	      if [[ $sline == *${suffix2}* ]];
                  then
                  sline=${sline%$suffix2}
              fi
              if [[ $sline == *${suffixss}* ]];
                  then
                  sline=${sline%$suffixss}
              fi
	      if [[ $sline == *${suffixhn}* ]];
                  then
                  sline=${sline%$suffixhn}
              fi
	      if [[ $sline == *${suffixhnfake}* ]];
                  then
                  sline=${sline%$suffixhnfake}
              fi
	      if [[ $sline == *${suffixhnfatjet}* ]];
              then
                  sline=${sline%$suffixhnfatjet}
              fi

	      if [[ ! -d "${sline2}" ]]; then
		  UNPROCESSED+=(${sline})
	      elif test "$(ls -A "$sline2")"; then
		  echo ${sline}
		  LISTOFSAMPLES+=(${sline})
	      else
		  UNPROCESSED+=(${sline})
	      fi
	  fi
	  
	  if [[ $submit_searchlist != "" ]];
	      then
	      if [[ $sline == *${submit_searchlist}* ]];
		  then 
		  prefix="SK"
		  suffix="_dilep"
		  suffixhn="_hndilep"
		  suffixhnfake="_hnfake"
		  suffixhnfatjet="_hnfatjet"
		  suffix2="_trilep"
		  suffixss="_sslep"
		  if [[ $sline == *${prefix}* ]];
		      then
		      sline=${sline:2}
		  fi
		  if [[ $sline == *${suffix}* ]];
		      then
		      sline=${sline%$suffix}
		  fi
		  if [[ $sline == *${suffix2}* ]];
                      then
                      sline=${sline%$suffix2}
                  fi
		  if [[ $sline == *${suffixss}* ]];
                      then
                      sline=${sline%$suffixss}
                  fi
		  if [[ $sline == *${suffixhn}* ]];
                  then
                      sline=${sline%$suffixhn}
		  fi
		  if [[ $sline == *${suffixhnfake}* ]];
                  then
                      sline=${sline%$suffixhnfake}
                  fi
		  if [[ $sline == *${suffixhnfatjet}* ]];
		  then
                      sline=${sline%$suffixhnfatjet}
                  fi

		  if [[ ! -d "${sline2}" ]]; then
		      UNPROCESSED+=(${sline})
		  elif test "$(ls -A "$sline2")"; then
		      echo ${sline}
		      LISTOFSAMPLES+=(${sline})
		  else
		      UNPROCESSED+=(${sline})
		  fi
	      fi
	  fi
      fi
      
    done < ${TXTPATH}"SNU_mc_"${submit_snuvlist}".txt"
    
    echo ""
    
    missing_comment="SNUanalyzer::sktree :: HELP :: If the sample you are looking for is not in the list above run 'sktree -A "$submit_snuvlist"'"$'\n'
    missing_comment+="SNUanalyzer::sktree :: HELP :: If 'sktree -A' shows sample is missing then we need to wait for the miniAOD to be produced\n"
    missing_comment+="SNUanalyzer::sktree :: HELP :: If 'sktree -A' shows sample is available at kisti then run 'sktree -r DATASETNAME' to request this sample\n"
    missing_comment+="SNUanalyzer::sktree :: HELP :: If 'sktree -A' shows sample is not there then no snuuple exists: run 'sktree -rsnu DATASETNAME' to request this snuuple"

    if [[ $isNoCut  == "true" ]];
        then

	counter=${#UNPROCESSED[@]}
        if [[ $counter -ne 0 ]];
            then
	    echo "Samples that have local flat snuuples but no NoCutskim are:"
	else  echo -e  $missing_comment
	fi
	for il in  ${UNPROCESSED[@]};
	  do
	  echo samplename = $il
	done
	echo ""
	if [[ $counter -ne 0 ]];
	    then
	    echo "If you want any of these sktrees run 'sktree -a SKTreeMakerNoCut -i <samplename> -c "$submit_snuvlist"'"
	fi
    fi
    if [[ $isLepton  == "true" ]];
        then
	counter=${#UNPROCESSED[@]}
	    if [[ $counter -ne 0 ]];
		then
		echo "Samples that have local flat snuuples but no lepton skim are:"
	    else    echo -e "$missing_comment"
	    fi
        for il in  ${UNPROCESSED[@]};
          do
          echo samplename = $il
        done
	echo ""
	if [[ $counter -ne 0 ]];
	    then
	    echo "If you want this sktree run 'sktree -a SKTreeMaker -i <samplename>  -c "$submit_snuvlist"'"
	    echo ""
	fi
    fi

    if [[ $isDiLep  == "true" ]];
        then
	counter=${#UNPROCESSED[@]}
	    if [[ $counter -ne 0 ]];
		then
		echo "Samples that have local flat snuuples but no dilepton skim are:"
	    else   echo -e $missing_comment

	    fi
	    for il in  ${UNPROCESSED[@]};
	      do
	      echo samplename = $il
	      
	    done
	    echo ""
	    if [[ $counter -ne 0 ]];
		then
		echo "If you want this sktree run 'sktree -a SKTreeMakerDiLep -i <samplename> -c "$submit_snuvlist"'"
		echo ""
	    fi
    fi

    
    if [[ $isTriLep  == "true" ]];
        then
        counter=${#UNPROCESSED[@]}
            if [[ $counter -ne 0 ]];
                then
                echo "Samples that have local flat snuuples but no trilepton skim are:"
            else   echo -e $missing_comment

            fi
            for il in  ${UNPROCESSED[@]};
              do
              echo samplename = $il

            done
            echo ""
            if [[ $counter -ne 0 ]];
                then
                echo "If you want this sktree run 'sktree -a SKTreeMakerTriLep -i <samplename> -c "$submit_snuvlist"'"
                echo ""
            fi
    fi
    
    if [[ $isSSLep  == "true" ]];
        then
        counter=${#UNPROCESSED[@]}
            if [[ $counter -ne 0 ]];
                then
                echo "Samples that have local flat snuuples but no sslepton skim are:"
            else   echo -e $missing_comment

            fi
            for il in  ${UNPROCESSED[@]};
              do
              echo samplename = $il

            done
            echo ""
            if [[ $counter -ne 0 ]];
		then
                echo "If you want this sktree run 'sktree -a SKTreeMakerssLep -i <samplename> -c "$submit_snuvlist"'"
		echo ""
            fi
    fi

    
    if [[ $specified_snuversion == "false" ]];
	then
	#Get number of snuversions
	for ic in  ${list_of_snuversions[@]};
	  do

	  if [[ $ic == $SNUVERSION ]];
	      then continue;
	  fi
	  if [[ $ic == *"v7-4"* ]];
	      then
	      echo "For snuversion "  $ic " the naming changed. Run 'sktree -L " $submit_skim " " $ic "'"
	      continue;
          fi

	  echo "Following samples are available in Snuversion: " ${ic} " but not in: " $SNUVERSION

	  if [[ $submit_skim  == "SKTree_NoSkim" ]];
	      then
	      check_path=$SKTREE_MC${ic}"/SKTrees/MCNoCut/"
	  fi
	  if [[ $submit_skim  == "SKTree_LeptonSkim" ]];
	      then
              check_path=$SKTREE_MC${ic}"/SKTrees/MC/"
          fi
	  if [[ $submit_skim  == "SKTree_DiLepSkim" ]];
	      then
	      check_path=$SKTREE_MC${ic}"/SKTrees/MCDiLep"
	  fi
	  if [[ $submit_skim  == "SKTree_TriLepSkim" ]];
              then
              check_path=$SKTREE_MC${ic}"/SKTrees/MCTriLep"
          fi
	  if [[ $submit_skim  == "SKTree_SSLepSkim" ]];
              then
              check_path=$SKTREE_MC${ic}"/SKTrees/MCSS"
          fi

	  while read line
	    do
	    if [[ $line == *${check_path}* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')

		if [[ $submit_searchlist == "" ]];
		    then
		    isDuplisnue=false
		    for il in  ${LISTOFSAMPLES[@]};
		      do

		      if [[ $sline == *${il}* ]];
			  then
			  isDuplisnue=true
		      fi
		    done
		    if [[ $isDuplisnue == "false" ]];
			then
			prefix="SK"
			suffix="_dilep"
			suffixhn="_hndilep"
			suffixhnfake="_hnfake"
			suffixhnfatjet="_hnfatjet"
			suffix2="_trilep"
			suffixss="_sslep"
			if [[ $sline == *${prefix}* ]];
			    then
			    sline=${sline:2}
			fi
			if [[ $sline == *${suffix}* ]];
			    then
			    sline=${sline%$suffix}
			fi
			if [[ $sline == *${suffix2}* ]];
                            then
                            sline=${sline%$suffix2}
                        fi
			if [[ $sline == *${suffixss}* ]];
			then
                            sline=${sline%$suffixss}
                        fi

			if [[ $sline == *${suffixhn}* ]];
			then
			    sline=${sline%$suffixhn}
			fi
			if [[ $sline == *${suffixhnfake}* ]];
                        then
                            sline=${sline%$suffixhnfake}
                        fi
			if [[ $sline == *${suffixhnfatjet}* ]];
                        then
			    sline=${sline%$suffixhnfatjet}
                        fi

			if [[ -d "${sline2}" ]]; then
			    if test "$(ls -A "$sline2")"; then
				echo ${sline}
				LISTOFSAMPLES+=(${sline})
			    fi
			fi
		    fi
		fi
		if [[ $submit_searchlist != "" ]];
		    then
		    isDuplisnue=false
		    if [[ $sline == *${submit_searchlist}* ]];
			then
			for il in  ${LISTOFSAMPLES[@]};
			  do
			  if [[ $sline == *${il}* ]];
			      then
			      isDuplisnue=true
			  fi
			done
			if [[ $isDuplisnue == "false" ]];
			    then
			    prefix="SK"
			    suffix="_dilep"
			    suffixhn="_hndilep"
			    suffixhnfake="_hnfake"
			    suffixhnfatjet="_hnfatjet"
			    suffix2="_trilep"
			    suffixss="_sslep"
			    if [[ $sline == *${prefix}* ]];
				then
				sline=${sline:2}
			    fi
			    if [[ $sline == *${suffix}* ]];
				then
				sline=${sline%$suffix}
			    fi
			    if [[ $sline == *${suffix2}* ]];
                                then
                                sline=${sline%$suffix2}
                            fi
			    if [[ $sline == *${suffixss}* ]];
							then
                                sline=${sline%$suffixss}
                            fi

			    if [[ $sline == *${suffixhn}* ]];
			    then
				sline=${sline%$suffixhn}
			    fi
			    if [[ $sline == *${suffixhnfake}* ]];
                            then
                                sline=${sline%$suffixhnfake}
                            fi
			    if [[ $sline == *${suffixhnfatjet}* ]];
			    then
                                sline=${sline%$suffixhnfatjet}
                            fi

			    if [[  -d "${sline2}" ]]; then
				
				if test "$(ls -A "$sline2")"; then
				    echo ${sline}
				    LISTOFSAMPLES+=(${sline})
				fi
			    fi
			fi
		    fi
		fi
	    fi
	  done < ${TXTPATH}"SNU_mc_"${ic}".txt"

	done
      fi

}				
    

while [ "$1" == "" ]; do
    usage
    exit
done

while [ "$1" != "" ]; do
    
    case $1 in
        -a | --analysis_name )  shift
                                submit_analyzer_name=$1
				set_submit_analyzer_name=true
				;;
        -i | --input )          shift
                                submit_file_tag=$1
				set_submit_file_tag=true
                                ;;
	-list | --inputlist )   shift
                                submit_file_list=$1
				set_submit_file_list=true
                                ;;
	
	-d | --debug_mode )     shift
	                        job_loglevel="$1"
				changed_job_loglevel=true
                                ;;
        -p | --data_period )    shift
                                job_data_lumi="$1"
				changed_job_data_lumi=true
                                ;;
        -F | --submitall)       shift
                                job_submitallfiles="true"
				;;
        -SIG | --submitall)     shift
                                RUNSIG="true"
				TXTPATH=${ANALYZER_RUN_PATH}"/txt/datasets_snu_sig_"
				;;
#	-sktree | --usesktrees )shift
#                                submit_skinput="$1"
#				changed_skinput=true
#                                ;;
        -s | --useskim )        shift
                                job_skim="$1"
				changed_skim=true
                                ;;
	-n | --njobs )          shift
	                        job_njobs=$1
				changed_job_njobs=true
				;;
        -ns| --njobs_skmaker )  shift
	                        set_sktreemaker_debug=true				
                                ;;

        -r | --requestSKtree )  shift
	                        request_sample=$1
	                        sendrequest
				exit 1
				;;
        -rsnu| --requestSNU )   shift
                                request_sample=$1
                                sendrequestsnu
                                exit 1
                                ;;        
	-S | --SampleTag  )     shift
                                submit_sampletag=$1
				set_submit_sampletag=true
                                ;;
	-q | --queue  )         shift
                                queuename=$1
                                ;;
        -qlist )                shift
                                listqueue
	                        exit 1
                                ;;

	-c | --SnuVersion)      shift
				submit_version_tag="$1"
				changed_submit_version_tag=true
				;;
        -ac | --AllSnuVersion)  shift
                                check_all_snuversions=$1
                                ;;
	-l | --file_tag_list)   shift
				submit_skim="FLATSNU"
				submit_searchlist=$1
				submit_snuvlist=$2
	                        runlist
				exit 1
				;;
        -L | --sktree_tag_list) shift
                                submit_skim=$1
                                submit_searchlist=$2
				submit_snuvlist=$3
                                runlist
                                exit 1
                                ;;
        -attachhist)                   shift
                                submit_draw=$1
                                ;;

        -A | --AvailableSnuuples) shift
	                        submit_snuvlist=$1
				submit_searchlist=$2
				listavailable
                                exit 1
                                ;;
        -t)                     shift
                                submit_dir_tag=$1
                                ;;

        -tagdiff)		shift
	                        submit_snu_tag=$1
				submit_snu_tag2=$2
				print_tag_diff
				exit 1
				;;
        -tagdiff7)              shift
                                submit_snu_tag=$1
                                submit_snu_tag2=$2
                                print_tag_diff7
                                exit 1
                                ;;

        -sktreelog)             shift
				submit_analyzer_name=$1
				submit_snuvlist=$2
				print_sktreemaker_logfile
				exit 1
                                ;;

	-updateselection)       shift
	                        object=$1
				update_selection
				exit 1
				;;
	-printID )              shift
	                        idname=$1
				idname2=$2
				printid
				exit 1
				;;
	-filename )             shift
	                        job_tmp_filename=$1
				;;
	-D | --GetProductionInfo)  shift
                                submit_snuvlist=$1
                                getalldatasetinfo
                                exit 1
                                ;;
        -userflag          )    shift
                                submit_skflag=$1
				;;


	-miniaod | --GetDataSetName)  shift
                            	submit_file_tag=$1
				submit_snuvlist=$2
				getdatasetname
				exit 1
				;;
        -M                      )  shift
                                make_sktrees=$1
                                ;;
	-V                      )  shift
	                           run_validation=$1
				   ;;

        -xsec | --GetDataSetXsec)  shift
                                submit_file_tag=$1
                                submit_snuvlist=$2
                                getdatasetxsec
                                exit 1
				;;
        -efflumi | --GetDataSetLumi)  shift
                                submit_file_tag=$1
                                submit_snuvlist=$2
                                getdatasetefflumi
                                exit 1
				;;
        -o | --output)          shift
	                        job_output_dir=$1
				changed_job_output_dir=true
                                ;;

        -G | --getoutputdir)       shift
                                GetOutPutDir="True"
                                submit_analyzer_name=$1
				submit_file_tag="WW"
				set_submit_file_tag=true
                                ;;

	-g | --file_tag_groups) shift
                                submit_searchlist=$1
                                rungroupedlist
				exit 1
                                ;;

        -tau | --run_tau_analyzer)    shift
                                job_run_taus=$1
                                 ;;
        -fake | --run_fake_analyzer)    shift
                                job_run_fake=$1
                                 ;;
        -m                  )    shift
                                submit_sk_message=$1
                                ;;
        -flip | --run_flip_analyzer) shift
                                job_run_flip=$1
                                ;;
        -b | --bkg ) shift
                                job_run_bkg=$1
                                ;;

        -events | --number_of_events) shift
	                        job_nevents=$1
				;;
	-nskip | --number_of_events) shift
                                job_nskip=$1
				;;

	-h | --help )        	shift
	                        usage
                         	other_commands=$1 
				if [[ $other_commands == "" ]];
				    then
				    echo "###########################Running command##################################################################################"
				    echo "Tag    |   Options                                            | DEFAULT PARAMETER   | COMMENT                             | "
				    echo "__________________________________________________________________________________________________________________________|"
				    echo "       |                                                      |                     |                                     |"
				    echo "-a     |   HNDiElectron/ExampleDiMuon/etc                     | default = ''        | Name of analysis class              |"
				    echo "       |   'sktree -a' lists options                          |                     |                                     |"
				    echo "       |                                                      |                     |                                     |"
				    echo "-S     |   ALL/DATA/MC/DoubleEG/DoubleMuon/                   | default = ''        | 'DATA' runs every data dataset.     |"
				    echo "       |   MuonEG/SingeMuon/SinglePhoton/SingleElectron       |                     | 'MC' runs every MC sample           |"
				    echo "       |                                                      |                     |                                     |"
				    echo "-s     |SKTree_NoSkim/SKTree_LeptonSkim/SKTree_Di[Tri]LepSkim | default='Lepton'    | Sets skim to use:                   | "
				    echo "       |  FLATSNU sets input to flatsnuuple not sktee         |                     | NoCuts/Lepton/DiLeptonstill work    | "
				    echo "       |                                                      |                     |                                     |"
				    echo "-n     |   #number of subjobs  (any number < 15)              | default=15          | default is 5 is -sktree=False       | "
				    echo "       |                                                      |                     |                                     |"
				    echo "-i     |   (i.e., DY10to50_MSnuNLO) run 'sktree -L' for more  | default = ''        | For running single MC samples:      | "
				    echo "       |                                                      |                     |                                     |"
				    echo "-list  |   (i.e., diboson_pythia) run 'sktree -g' for more    | default = ''        | For running on list of MC samples.  | "
				    echo "       |                                                      |                     |                                     |"
				    echo "-c     |   snuversion of inputfile                            | default = ${SNUVERSION}    | (only needed if not running         |"
				    echo "       |                                                      |                     | default/latest)                     | " 
                                    echo "-ac    |   true/false                                         | default = 'false'   | Check all snuversions for input     | "
                                    echo "       |                                                      |                     | Rare: only set true if a sample is  |" 
                                    echo "       |                                                      |                     | not available in  ${SNUVERSION}            |"
                                    echo "       |                                                      |                     |                                     |" 
				    echo "-d     |   debug mode : INFO/DEBUG/WARNING                    | default = INFO      |                                     | "
				    echo "       |                                                      |                     |                                     |"
				    echo "-p     |   period to run in data/normalise inMC: C/D/CtoD(ALL)| default = CtoD      | Only change if you wish to run      | "
				    echo "       |                                                      |                     | on a single data period             | "
                                    echo "-o     |   setoutput directory.                               | default= ''         | Does not work for SKTreeMaker Code. |"
                                    echo "       |                                                      |                     |                                     |"
                                    echo "-fake  |   set flag for fake analyzer                         | default = 'false'   |                                     |"
                                    echo "       |                                                      |                     |                                     |"
                                    echo "-flip  |   set flag for charge flip analyzer                  | default = 'false'   |                                     |"
                                    echo ""
				    echo "Run 'sktree -h more' or 'sktree -h  debug' for more commands"
				fi
				
				if [[ $other_commands == "debug" ]];
				    then
                                    echo "###########################Other command#####################################################################################"
                                    echo "Tag    |   Options                                            | DEFAULT PARAMETER   | COMMENT                             | "
                                    echo "__________________________________________________________________________________________________________________________|"
				    echo "-events|   number of events to process                        |  default = ''       |                                     |"
				    echo "-nskip |   number of events to skip                           |  default = ''       |                                     |"  
				    echo "       |                                                      |                     |                                     |"
				    
				fi
				
				if [[ $other_commands == "more" ]];
	  			    then
				    echo "###########################Other command#####################################################################################"
				    echo "Tag      |   Options                                          | DEFAULT PARAMETER   | COMMENT                             | "
				    echo "__________________________________________________________________________________________________________________________|"
				    echo "-l       | (can give search/snuversion as an option )         | default = ''        | returns a list of available         | "
				    echo "         | i.e.  sktree -l QCD  OR  sktree -l DY v7-6-3       |                     | datasets in each snuversion         |"
				    echo "-L       | can give search/skim/snuversion as on option)      | default = ''        | returns a list of available         | "
				    echo "         | use like 'sktree -L DiLep QCD                      |                     | sktrees  in each snuversion         |" 
				    echo "         |                                                    |                     |                                     |"
				    echo "-g       |                                                    | default = ''        | returns a list of available input   |"
				    echo "         |                                                    |                     | arrays to input with:               |"
				    echo "         |                                                    |                     | sktree -list command                |" 
				    echo "         |                                                    |                     | Only available from v7-6-3          |"
                                    
				    echo "-tagdiff | tagname                                            | default = ''        | print change log between current tag|"
				    echo "         |                                                    |                     | and <tagname>                       |"
				    echo "-sktreelog| class name : vatversion                           | default = ''        | print log of sktreemaker            |"
				    echo "         |                                                    |                     | Only available from v7-6-4          |"

				    echo "-D       | any allowed  SNUVERSION                            | default = $SNUVERSION    | returns Info on Snuuple production.|"
				    echo "-G       | any allowed  Analyzer                              | default = ''        | returns drfault outputdir.|"
				    echo "-printID | any allowed  ID name  (i.e., MUON_POG_TIGHT)       | default = ''        | returns Info on object id.|"
				    echo "-userflag| Get user flag   flag1,flag2                        | default = ''        |  pass in string                     |"
				    
				    
				    echo "-miniaod | file_tag (i.e., DY10to50_MSnuNLO)                  | default = ''        | returns datasetname.                |"
				    echo "         |                                                    |                     | Only available from v7-6-3          |"  
                                    echo "-xsec    | file_tag (i.e., DY10to50_MSnuNLO)                  | default = ''        | returns dataset xsec                |"
                                    echo "         |                                                    |                     | Only available from v7-6-3          |"
                                    echo "-efflumi | file_tag (i.e., DY10to50_MSnuNLO)                  | default = ''        | returns dataset lumi.               |"
                                    echo "         |                                                    |                     | Only available from v7-6-3          |"

				    echo "-r       | sktree -A to see possible samples                  | default = ''        | sends email request to make sktree  |"
				    echo "         |                                                    |                     | that is not current available at snu|"
				    echo "         |                                                    |                     |                                     |"
				    echo "-updateselection | any allowed objectname                           | default = ""  | updates selection file and sends email.|"
				    echo "-A       | can speficy snuversion and search                  | default = ${SNUVERSION}    | lists missing samples due to no     |"
				    echo "         |                                                    |                     | MiniAOD and available samples       |"
				    echo "         |                                                    |                     | that can be processed               |"        
				    echo " -attachhist        | set True or False                                            |                     | that can be processed               |"        
				    
				fi
				exit
                                ;;
        * )                     usage
	exit 1
    esac
    shift
done


############################################################



#declare -a streams=("")
declare -a data_periods=("")
declare -a ALL=("DoubleMuon" "DoubleEG" "MuonEG" "SinglePhoton" "SingleElectron" "SingleMuon" "DoubleMuon_CF")


if [[ $job_data_lumi == "ALL" ]];
    then

    if [[ $SNUVERSION == "v8-0-7" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-6" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi

    if [[ $SNUVERSION == "v8-0-4" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-3" ]];then
	declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-2" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G")
    fi

    if [[ $SNUVERSION == "v8-0-1" ]];then
        declare -a data_periods=("B" "C" "D" "E")
    fi
    if [[ $SNUVERSION == "v7-6-6" ]];then
        declare -a data_periods=("C" "D")
    fi

fi

if [[ $job_data_lumi == $snudatatag  ]];
then
    if [[ $SNUVERSION == "v8-0-7" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-6" ]];then
         declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-4" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-3" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G" "H_v2" "H_v3")
    fi
    if [[ $SNUVERSION == "v8-0-2" ]];then
        declare -a data_periods=("B" "C" "D" "E" "F" "G")
    fi
    if [[ $SNUVERSION == "v8-0-1" ]];then
        declare -a data_periods=("B" "C" "D" "E")
    fi
    if [[ $SNUVERSION == "v7-6-6" ]];then
        declare -a data_periods=("C" "D")
    fi
fi

if [[ $job_data_lumi == "CtoD" ]];
    then
    declare -a data_periods=("C" "D")
    export snudatatag="CtoD"
fi
if [[ $job_data_lumi == "BtoE" ]];
    then
    declare -a data_periods=("B" "C" "D" "E")
    export snudatatag="BtoE"
fi


export SNUAnalyzerPeriod="None"

if [[ $job_data_lumi == "B" ]];
    then
    declare -a data_periods=("B")
    export SNUAnalyzerPeriod="B"
    export snudatatag="B"
	
fi

if [[ $job_data_lumi == "C" ]];
    then
    declare -a data_periods=("C")
    export SNUAnalyzerPeriod="C"
    export snudatatag="C"
fi
if [[ $job_data_lumi == "D" ]];
    then
    declare -a data_periods=("D")   
    export SNUAnalyzerPeriod="D"                                                                                        
    export snudatatag="D"

fi 
if [[ $job_data_lumi == "E" ]];
    then
    declare -a data_periods=("E")
    export SNUAnalyzerPeriod="E"                                                                                        
    export snudatatag="E"
    
fi
if [[ $job_data_lumi == "F" ]];
    then
    declare -a data_periods=("F")
    export SNUAnalyzerPeriod="F"                                                                                        
    export snudatatag="F"

fi
if [[ $job_data_lumi == "G" ]];
    then
    declare -a data_periods=("G")
    export SNUAnalyzerPeriod="G"
    export snudatatag="G"
	
fi
if [[ $job_data_lumi == "H" ]];
    then
    declare -a data_periods=("H_v2" "H_v3")
    export SNUAnalyzerPeriod="H"
    export snudatatag="H"
	
fi
if [[ $job_data_lumi == "H_v3" ]];
    then
    declare -a data_periods=( "H_v3")
    export SNUAnalyzerPeriod="H"
    export snudatatag="H"

fi

if [[ $job_data_lumi == "H_v2" ]];
    then
    declare -a data_periods=("H_v2")
    export SNUAnalyzerPeriod="H"
    export snudatatag="H"
fi

if [[ $job_data_lumi == "GH" ]];
    then
    declare -a data_periods=("G" "H_v2" "H_v3")
    export SNUAnalyzerPeriod="GH"
    export snudatatag="GH"
	
fi

ARG1=snudataperiods
eval getlist_cv=(\${$ARG1[@]})
for dataperiod in  ${getlist_cv[@]};
do
    if [[ $job_data_lumi == $dataperiod ]];
    then
    declare -a data_periods=($dataperiod)
    export SNUAnalyzerPeriod=$dataperiod
    fi
done


if [[ $submit_file_tag  != "" ]];
    then
    runMC=true
fi


if [[ $submit_file_list  != "" ]];
    then
    runMC=true
fi


if [[ $submit_sampletag  == "ALL" ]];
    then
    runMC=true
    submit_file_list="all_mc"
    runDATA=true
fi
if [[ $submit_sampletag  == "MC" ]];
    then
    runMC=true
    submit_file_list="all_mc"

fi

declare -a  DATA=("DoubleMuon" "DoubleEG" "MuonEG"  "SingleElectron" "SingleMuon")
if [[ $submit_sampletag  == "DATA" ]];
    then
    runDATA=true


fi    
declare -a  DATADILEP=("DoubleMuon" "DoubleEG" "MuonEG")
if [[ $submit_sampletag  == "DATADILEP" ]];
    then
    runDATA=true
fi



declare -a DoubleEG=("DoubleEG")

if [[ $submit_sampletag  == "DoubleEG" ]];
    then
    runDATA=true
fi

declare -a DoubleMuon=("DoubleMuon")
if [[ $submit_sampletag  == "DoubleMuon" ]];
    then
    runDATA=true
fi

declare -a DoubleMuon_CF=("DoubleMuon_CF")
if [[ $submit_sampletag  == "DoubleMuon_CF" ]];
    then
    runDATA=true
fi


declare -a MuonEG=("MuonEG")
if [[ $submit_sampletag  == "MuonEG" ]];
    then
    runDATA=true
    
fi

declare -a SinglePhoton=("SinglePhoton")
if [[ $submit_sampletag  == "SinglePhoton" ]];
    then
    runDATA=true
fi

declare -a SingleElectron=("SingleElectron")

if [[ $submit_sampletag  == "SingleElectron" ]];
    then
    runDATA=true
fi

declare -a SingleLepton=("SingleElectron" "SingleMuon")

if [[ $submit_sampletag  == "SingleLepton" ]];
    then
    runDATA=true
fi



declare -a SingleMuon=("SingleMuon")

if [[ $submit_sampletag  == "SingleMuon" ]];
    then
    runDATA=true
fi

if [[ $submit_version_tag == *"v7-4"* ]];
    then
    declare -a SingleElectron=("singleelectron")
    declare -a SingleMuon=("singlemuon")
    declare -a DoubleMuon=("muon")
    declare -a MuonEG=("emu")
    declare -a DoubleEG=("egamma")
fi

#if [[ $submit_version_tag != "" ]];
#    then
#    #checkdata
#fi


declare -a FULLLISTOFSAMPLES=()
declare -a FULLLISTOFSAMPLESNOCUT=()
declare -a FULLLISTOFSAMPLESLEPTON=()
declare -a FULLLISTOFSAMPLESDILEP=()
declare -a FULLLISTOFSAMPLESTRILEP=()
declare -a FULLLISTOFSAMPLESTSSLEP=()


if [[ $check_all_snuversions != "true" ]];
    then
    if [[ $check_all_snuversions != "false" ]];
	then
	if [[ $check_all_snuversions == "" ]];
	    then
	    echo "No input for -ac command: this should be '-ac true', since false is default"
	    exit 1
	fi
	
	echo "wrong setting for -ac command: this should be '-ac true', since false is default" 
	exit 1
    fi
fi
 
####  MAKING FULL LIST OF SAMPLES. THIS IS TO CHECK THE INPUT IS CORRECT
### Therefore we only make the lists if there is some input
MakeFullLists=false
if [[ $submit_file_list  != ""  ]];
    then
    MakeFullLists=true
fi
if [[ $submit_file_tag  != ""  ]];
    then
    MakeFullLists=true
fi  
### We make only the list for the skim specified: To check the sample is available in with this skim
### INCASE skim is set wrong we make a quick check
if [[ $changed_skim != "true" ]];
    then
    if [[ $changed_skinput  == "true" ]];
	then
        if [[ $submit_skinput == "false" ]];
            then
            job_skim=FLATSNU
        fi
    fi
fi

#### Since SKTreeMaker codes have skims set as default: Never need to change we set this here
if [[ $submit_analyzer_name == *"SKTreeMaker"* ]];
    then
    job_skim=FLATSNU
    if [[ $submit_analyzer_name == "SKTreeMakerDiLep" ]];
	then
	job_skim=SKTree_LeptonSkim
    fi
    if [[ $submit_analyzer_name == "SKTreeMakerTriLep" ]];
        then
        job_skim=SKTree_LeptonSkim
    fi
    if [[ $submit_analyzer_name == "SKTreeMakerSSLep" ]];
        then
	    job_skim=SKTree_LeptonSkim
    fi


fi


####### If full list is needed
if [[ $MakeFullLists == "true" ]];
    then
    
    ##### LOOP OVER ALL SNUVERSIONS (ONLY STORE SAMPLES THAT EXIST AND NOT IN NEWER VERSIONS)
    for iclist in  ${list_of_snuversions[@]};
      do
      #### IF CHECK ALL SNUVERSION  = FALSE WE LOOK AT ONLY $SNUVERSION
      if [[ $check_all_snuversions != "true" ]];
	  then
	  if [[ $iclist != ${submit_version_tag} ]];
	      then
	      continue
	  fi
      fi

      
      #### LOOP OVER INPUT TXT FILE AND CHECK FOR AVAILABLE SAMPLES
      
      while read line
	do
	
	if [[ $job_skim == "FLATSNU" ]];
	    then
	    
	    if [[ $line == *$FLATSNU_MC* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')
		
		isDuplisnue=false
		for il in  ${FULLLISTOFSAMPLES[@]};
		  do
		  if [[ $sline == $il ]];
		      then
		      isDuplisnue=true
		  fi
		done
		if [[ $isDuplisnue == "false" ]];
		    then
		    
		    if [[ -d "${sline2}" ]]; then
			if test "$(ls -A "$sline2")"; then
			    FULLLISTOFSAMPLES+=(${sline})
			fi
		    fi
		fi
	    fi
	fi
	
	if [[ $job_skim == "SKTree_NoSkim" ]];
	    then
	    if [[ $line == *$SKTREE_MC${iclist}"/SKTrees/MCNoCut/"* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')
		
		prefix="SK"
		suffix="_nocut"
		if [[ $sline == *${prefix}* ]];
		    then
		    sline=${sline:2}
		fi
		if [[ $sline == *${suffix}* ]];
		    then
		    sline=${sline%$suffix}
		fi
		
		
		isDuplisnue=false
		for il in  ${FULLLISTOFSAMPLESNOCUT[@]};
		  do
		  
		  if [[ $sline == $il ]];
		      then
		      isDuplisnue=true
		  fi
		done
		if [[ $isDuplisnue == "false" ]];
		    then
		    if [[ -d "${sline2}" ]]; then
			if test "$(ls -A "$sline2")"; then
			    FULLLISTOFSAMPLESNOCUT+=(${sline})
			fi
		  fi
		fi
	    fi
	fi

	if [[ $job_skim == "SKTree_LeptonSkim" ]];
	    then	
	    
	    checkline=$SKTREE_MC${iclist}"/SKTrees/MC/"
	    if [[ ${iclist} == *"v7-4-4"* ]];
                then
                checkline=$SKTREE_MC${iclist}"/SKTrees/Sep15/MC/"
            fi

            if [[ $line == *$checkline* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')
		
		prefix="SK"
		if [[ $sline == *${prefix}* ]];
		    then
		    sline=${sline:2}
		fi
		
		isDuplisnue=false
		for il in  ${FULLLISTOFSAMPLESLEPTON[@]};
		  do
		  if [[ $sline == $il ]];
		      then
		      isDuplisnue=true
		  fi
		done
		if [[ $isDuplisnue == "false" ]];
		    then
		    if [[ -d "${sline2}" ]]; then
			if test "$(ls -A "$sline2")"; then
			    FULLLISTOFSAMPLESLEPTON+=(${sline})
			fi
		    fi
		fi
	    fi
	fi  
	if [[ $job_skim == "SKTree_DiLepSkim" ]];
	    then
	    checkline=$SKTREE_MC${iclist}"/SKTrees/MCDiLep"
	    if [[ ${iclist} == *"v7-4-4"* ]];
		then
		checkline=$SKTREE_MC${iclist}"/SKTrees/Sep15/MCDiLep"
	    fi


	    if [[ $line == *$checkline* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')
		
		prefix="SK"
		suffix="_dilep"
		if [[ $sline == *${prefix}* ]];
		    then
		    sline=${sline:2}
		fi
		if [[ $sline == *${suffix}* ]];
		    then
		    sline=${sline%$suffix}
		fi
		
		isDuplisnue=false
		for il in  ${FULLLISTOFSAMPLESDILEP[@]};
		do
		  if [[ $sline == $il ]];
		      then
		      isDuplisnue=true
		  fi
		done
		if [[ $isDuplisnue == "false" ]];
		    then
		    if [[ -d "${sline2}" ]]; then
			if test "$(ls -A "$sline2")"; then
			    FULLLISTOFSAMPLESDILEP+=(${sline})
			fi
		    fi
		fi
	    fi
	fi    
	if [[ $job_skim == "SKTree_TriLepSkim" ]];
	    then
	    checkline=$SKTREE_MC${iclist}"/SKTrees/MCTriLep"
	    if [[ ${iclist} == *"v7-4-4"* ]];
		then
		checkline=$SKTREE_MC${iclist}"/SKTrees/Sep15/MCTriLep"
	    fi
	    
	    
	    if [[ $line == *$checkline* ]];
		then
		sline=$(echo $line | head -n1 | awk '{print $1}')
		sline2=$(echo $line | head -n1 | awk '{print $6}')
		
		prefix="SK"
		suffix="_trilep"
		if [[ $sline == *${prefix}* ]];
		    then
		    sline=${sline:2}
		fi
		if [[ $sline == *${suffix}* ]];
		    then
		    sline=${sline%$suffix}
		fi
		
		isDuplisnue=false
		for il in  ${FULLLISTOFSAMPLESTRILEP[@]};
		  do
		  if [[ $sline == $il ]];
			  then
		      isDuplisnue=true
		  fi
		done
		if [[ $isDuplisnue == "false" ]];
		    then
		    if [[ -d "${sline2}" ]]; then
			if test "$(ls -A "$sline2")"; then
			    FULLLISTOFSAMPLESTRILEP+=(${sline})
			fi
		    fi
		fi
	    fi
	fi
	if [[ $job_skim == "SKTree_SSLepSkim" ]];
        then
            checkline=$SKTREE_MC${iclist}"/SKTrees/MCSS"
            if [[ ${iclist} == *"v7-4-4"* ]];
                then
                checkline=$SKTREE_MC${iclist}"/SKTrees/Sep15/MCSS"
            fi


            if [[ $line == *$checkline* ]];
                then
                sline=$(echo $line | head -n1 | awk '{print $1}')
                sline2=$(echo $line | head -n1 | awk '{print $6}')

                prefix="SK"
                suffix="_SSlep"
                if [[ $sline == *${prefix}* ]];
                    then
                    sline=${sline:2}
                fi
                if [[ $sline == *${suffix}* ]];
                    then
                    sline=${sline%$suffix}
                fi

                isDuplisnue=false
                for il in  ${FULLLISTOFSAMPLESSS[@]};
                  do
                  if [[ $sline == $il ]];
                          then
                      isDuplisnue=true
                  fi
                done
                if [[ $isDuplisnue == "false" ]];
                    then
                    if [[ -d "${sline2}" ]]; then
                        if test "$(ls -A "$sline2")"; then
                            FULLLISTOFSAMPLESSS+=(${sline})
                        fi
                    fi
                fi
            fi
        fi

      done < ${TXTPATH}"SNU_mc_"${iclist}".txt"
    done
    
fi

