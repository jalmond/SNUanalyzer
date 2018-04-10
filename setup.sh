#! /bin/sh
# $Id: setup.sh 1 01/12/2013 jalmond $
###################################################################################
# @Project: Analyzer/SKTree - ROOT-based analysis framework for Korea CMS group #
#                                                                                 #
# @author John Almond       <jalmond@cern.ch>           - SNU                     #
#                                                                                 #
# Script that has to be sourced before compiling/running SFrame.                  #
#                                                                                 #
##################################################################################

# Greet the user

echo "Setting up environment for compiling/running CATAnalzer with SKTree"

setupok=False

Flag=""

if [[ $1 == "" ]]; then
    Flag="SNU"

else
    Flag="LQ"
fi

function killbkg {
    python $PWD/python/killbkg.py -i $1
}

while read line
  do
  if [[ $line == *"LANG"* ]]; then
      setupok=True
  fi
done < ~/.bash_profile

if [[ $setupok == "False" ]]; then
    echo "Please add the following lines to ~/.bash_profile file:"
    echo 'LANG="en_US.utf-8"'
    echo 'LC_COLLATE="lt_LT.utf-8"'
    echo 'LC_TIME="en_DK.utf-8"'
    return 1
fi

if [[ $USER == "jalmond" ]]; then
    alias cat_path_analysis_ls='ll -rth /data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/HNDiLepton/periodBtoH/ '
    if [ $LQANALYZER_DIR ]; then
	echo "Running on batch"
    else
	source python/jalmondsetup.sh
    fi
    function cat_path_analysis_ls {
        ll -rth  /data2/CAT_SKTreeOutput/JobOutPut/jalmond/${Flag}analyzer/data/output/CAT/HNDiLepton/periodBtoH/${1}
    }
    function cat_path_analysis {
	cd /data2/CAT_SKTreeOutput/JobOutPut/jalmond/${Flag}analyzer/data/output/CAT/HNDiLepton/periodBtoH/${1}
    }
fi




if [[ $PWD !=  *"/data4/${Flag}AnalyzerCode/"* ]];
then
    if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
    then
        echo "Setup failed. ${Flag}analyzer needs to be in /data4/${Flag}AnalyzerCode/"$USER
        if [ ! -d /data4/${Flag}AnalyzerCode/$USER ]; then
            mkdir /data4/${Flag}AnalyzerCode/$USER
        fi
        echo "Move the current ${Flag}Analyzer directory to "/data4/${Flag}AnalyzerCode/$USER

        return
    fi
fi


if [ $LQANALYZER_DIR ]; then
    echo ${Flag}ANALYZER_DIR is already defined, use a clean shell
    return 1
fi






## variables that are specific to your machine: Change if noy listed
if [ "$HOSTNAME" = "cms2.snu.ac.kr" ] || [ "$HOSTNAME" = "cms1.snu.ac.kr" ]; then    
    source /share/apps/root_v5-34-32/root/bin/thisroot.sh 
elif [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    source /share/apps/root_v5_34_32/root/bin/thisroot.sh
else
    source /share/apps/root_v5-34-32/root/bin/thisroot.sh
fi    


# speficy the ${Flag}ANALYZER_DIR base directory, i.e., the directory in which this file lives
export ANALYZER_DIR=${PWD}



if [[ $USER == "jalmond" ]]; then
    python python/setupAN.py
fi



if [[ $1 == *"v7"* ]]; then
    echo "Setting up tag "$1
    export CHECKTAGFILE=/data1/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/SetBrachAndTag_$1.sh
    if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
    then
	export CHECKTAGFILE=/data2/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/SetBrachAndTag_$1.sh
    fi

    if [[ ! -f $CHECKTAGFILE ]]; then 
	export ANALYZER_DIR=""
	echo $1 "is not allowed input. Use one of:"
	
	source /data1/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/$CATVERSION.sh
	if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
	then
	    source /data2/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/$CATVERSION.sh
	fi
	for ic in  ${list_of_catversions[@]};
        do
            echo $ic
	done
	return 1
    
    fi
    export ANALYZER_MOD="/data1/${Flag}Analyzer_rootfiles_for_analysis/CATMOD2015/"
    if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
        then
	export ANALYZER_MOD="/data2/${Flag}Analyzer_rootfiles_for_analysis/CATMOD2015/"

    fi
    source $ANALYZER_DIR/bin/setup2015.sh
    export running2015=True
    cvdir=$ANALYZER_DIR/${Flag}Lib/$CATVERSION
    if [[ ! -d "${cvdir}" ]]; then
        mkdir $cvdir
        make distclean
	make
    fi
    return 1
fi

export ANALYZER_MOD="/data1/${Flag}Analyzer_rootfiles_for_analysis/CATMOD/"
if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    export ANALYZER_MOD="/data2/${Flag}Analyzer_rootfiles_for_analysis/CATMOD/"

fi
python ${ANALYZER_DIR}/scripts/CheckEmailIsSetup.py
cat_email="NULL"
while read line
do
    prefix="email = "
    if [[ $line == $prefix* ]];
    then
        line=${line:${#prefix}}
        cat_email=$line
    fi
done < ${ANALYZER_DIR}/bin/catconfig
if [[ $cat_email  == "NULL" ]];
then
    echo "Email not setup. run setup.sh again"
    export ANALYZER_DIR=""
    return 1
fi

##### Check that this is not the branch and a tag was checked out
export CHECKTAGFILE=$ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
source $CHECKTAGFILE branch

source $ANALYZER_DIR/bin/CheckTag.sh

buglist=/data1/${Flag}Analyzer_rootfiles_for_analysis/CATTag/BuggyTag.txt
if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    buglist=/data2/${Flag}Analyzer_rootfiles_for_analysis/CATTag/BuggyTag.txt
fi

while read line
do
    if [[ $line == $CATTAG* ]];
    then
	echo "Current tag is buggy. Please update to newer tag."
        exit
    fi
done < $buglist

export LIBTAG=""
if [[ $1 != "" ]];then

    export CHECKTAGFILE=/data1/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/SetBrachAndTag_$1.sh
    if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
    then
	export CHECKTAGFILE=/data2/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/SetBrachAndTag_$1.sh
    fi
    if [[ ! -f $CHECKTAGFILE ]]; then
	export ANALYZER_DIR=""
        echo $1 "is not allowed input. Use one of:"
	
	if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
	then
            source /data2/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/$CATVERSION.sh
	else
	    source /data1/${Flag}Analyzer_rootfiles_for_analysis/CattupleConfig/$CATVERSION.sh

	fi
        for ic in  ${list_of_catversions[@]};
        do
            echo $ic
	done
        return 1
    fi
    
    source $CHECKTAGFILE branch
    export LIBTAG=$CATVERSION
fi

export yeartag="80X/"



alias cathistcounter="source scripts/Counter.sh "
alias catcutflowcounter="source scripts/CutFlow.sh "
alias sktree="bash submitSKTree.sh"
alias sktreemaker="bash submitSKTree.sh -M True "
alias sktree_val="bash submitSKTree.sh -V True "
alias sktree_bkg="nohup bash submitSKTree.sh -b True "
alias new_git_tag="bash "$ANALYZER_DIR"/scripts/setup/git_newtag.sh"
alias git_commit_lq="bash scripts/setup/git_commit.sh"
alias sktree_bkg_log="python python/PrintBkgJob.py"

export ANALYZER_FILE_DIR="/data1/${Flag}Analyzer_rootfiles_for_analysis/CATAnalysis2016/"
export ANALYZER_DATASETFILE_DIR="/data1/${Flag}Analyzer_rootfiles_for_analysis/DataSetLists/AnalysisFiles/"
export ANALYZER_DATASET_DIR="/data1/${Flag}Analyzer_rootfiles_for_analysis/DataSetLists/"
export ANALYZER_SKTreeLOG_DIR="/data1/${Flag}Analyzer_rootfiles_for_analysis/CATSKTreeMaker/"
export CATTAGDIR="/data1/${Flag}Analyzer_rootfiles_for_analysis/CATTag/"
if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    export ANALYZER_FILE_DIR="/data2/${Flag}Analyzer_rootfiles_for_analysis/CATAnalysis2016/"
    export ANALYZER_DATASETFILE_DIR="/data2/${Flag}Analyzer_rootfiles_for_analysis/DataSetLists/AnalysisFiles/"
    export ANALYZER_DATASET_DIR="/data2/${Flag}Analyzer_rootfiles_for_analysis/DataSetLists/"
    export ANALYZER_SKTreeLOG_DIR="/data2/${Flag}Analyzer_rootfiles_for_analysis/CATSKTreeMaker/"
    export CATTAGDIR="/data2/${Flag}Analyzer_rootfiles_for_analysis/CATTag/"
fi

export running2015=False

# Modify to describe your directory structure.
# all directories are below the ${Flag}Analyser base directory specified above
### setup paths to be used in analysis code
export ANALYZER_ANALYSIS_PATH=${ANALYZER_DIR}/${Flag}Analysis/
export ANALYZER_SRC_PATH=${ANALYZER_DIR}/${Flag}Analysis/Analyzers/src/
export ANALYZER_INCLUDE_PATH=${ANALYZER_DIR}/${Flag}Analysis/Analyzers/include/
export ANALYZER_CORE_PATH=${ANALYZER_DIR}/${Flag}Core/

export isSLC5="False"
export BTAGDIR=${ANALYZER_DIR}/${Flag}Analysis/AnalyzerTools/BTag/BTagC11/
export ROCHDIR=${ANALYZER_DIR}/${Flag}Analysis/AnalyzerTools/rochcor2016/
if [[ "$HOSTNAME" == "cms.snu.ac.kr" ]];
then 
    if [[ $LIBTAG == *"v"* ]]; then
	export OBJ=obj/cms21$LIBTAG
	export ANALYZER_LIB_PATH=${ANALYZER_DIR}/${Flag}Lib/cms21$LIBTAG/
        export ANALYZER_BATCHLIB_PATH=${ANALYZER_DIR}/${Flag}Lib/batch/
	
    else 
	export OBJ=obj/cms21
        export ANALYZER_LIB_PATH=${ANALYZER_DIR}/${Flag}Lib/cms21/
	export ANALYZER_BATCHLIB_PATH=${ANALYZER_DIR}/${Flag}Lib/batch/
    fi
elif [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    export OBJ=obj/cluster/
    export ANALYZER_LIB_PATH=${ANALYZER_DIR}/${Flag}Lib/cluster/

elif [[ "$HOSTNAME" == "cms1" ]];
then
    export OBJ=obj/cms1
    export ANALYZER_LIB_PATH=${ANALYZER_DIR}/${Flag}Lib/cms1/

else
    export OBJ=obj/cms2
    export ANALYZER_LIB_PATH=${ANALYZER_DIR}/${Flag}Lib/cms2/

fi

export ANALYZER_OLDLIB_PATH=${ANALYZER_DIR}/${Flag}Lib/

export ANALYZER_RUN_PATH=${ANALYZER_DIR}/${Flag}Run/
export ANALYZER_CLUSTER_TXT_PATH=${ANALYZER_DIR}/${Flag}Run/txt/Cluster/
export ANALYZER_BIN_PATH=${ANALYZER_DIR}/bin/
### set SKTree path
export SKTREE_INCLUDE_PATH=${ANALYZER_DIR}/${Flag}Core/SKTree/include/
## setup directory to store analysis rootfiles
export FILEDIR=${ANALYZER_DIR}/data/rootfiles/
export IDFILEDIR=${ANALYZER_DIR}/data/ID/80X/
export LUMIFILEDIR=${ANALYZER_DIR}/data/Luminosity/80X/
export TRIGGERFILEDIR=${ANALYZER_DIR}/data/Trigger/80X/
export BTAGFILEDIR=${ANALYZER_DIR}/data/BTag/80X/
export PILEUPFILEDIR=${ANALYZER_DIR}/data/Pileup/80X/



if [ ! -d ${ANALYZER_OLDLIB_PATH} ]; then
    echo Directory ${ANALYZER_OLDLIB_PATH} does not exist ... creating it
    mkdir ${ANALYZER_OLDLIB_PATH}
fi

if [ ! -d ${ANALYZER_LIB_PATH} ]; then
    echo Directory ${ANALYZER_LIB_PATH} does not exist ... creating it
    mkdir ${ANALYZER_LIB_PATH}
    file="${ANALYZER_OLDLIB_PATH}/libAnalysisCore.so"
    if [ -f "$file" ]; then
	echo Old lib dir ${ANALYZER_OLDLIB_PATH} is redundant. Will remove these library
	rm  ${ANALYZER_OLDLIB_PATH}/*.so
	rm  ${ANALYZER_OLDLIB_PATH}/*map
	rm  ${ANALYZER_CORE_PATH}/*/obj/*.o
	rm -r ${ANALYZER_CORE_PATH}/*/obj/dep/
	rm  ${ANALYZER_ANALYSIS_PATH}/*/obj/*.o
	rm -r ${ANALYZER_ANALYSIS_PATH}/*/obj/dep/
    fi
fi

### Load useful functions
source ${ANALYZER_BIN_PATH}/cleanup.sh 
### make directories that git does not allow to store

export ANALYZER_OUTPUT_PATH=/data2/CAT_SKTreeOutput/JobOutPut/${USER}/${Flag}analyzer/data/output/

export ANALYZER_LOG_PATH=/data2/CAT_SKTreeOutput/JobOutPut/${USER}/${Flag}analyzer/data/logfiles/
export ANALYZER_LOG_8TeV_PATH=${ANALYZER_DIR}/data/logfiles/

if [ $HOSTNAME == "tamsa2.snu.ac.kr" ];
then
    export ANALYZER_OUTPUT_PATH=/data4/CAT_SKTreeOutput/JobOutPut/${USER}/${Flag}analyzer/data/output/
    export ANALYZER_LOG_PATH=/data4/CAT_SKTreeOutput/JobOutPut/${USER}/${Flag}analyzer/data/logfiles/
fi

python ${ANALYZER_DIR}/python/SetUpWorkSpace.py
python ${ANALYZER_DIR}/python/BackUpDirectory.py
python ${ANALYZER_DIR}/python/SetupEmailList.py

# Setup root area and other paths
 
if [[ `which root-config` == "" ]]; then
    echo "Warning: ROOT environment doesn't seem to be configured!"
    source $root_setup
    if [[ `which root-config` == "" ]]; then
	echo  "Error: ROOT environment cannot be configured!"
    else echo "Setup root enviroment " 
    fi
fi

if [ -z ${ROOTSYS} ] ; then
    echo "Warning: ROOT environment doesn't seem to be configured!"
    echo "Add these lines to your ~/.bashrc file to remove this warning in future."
    echo ""
    echo "source /usr/local/bin/thisroot.sh"
    echo ""
    export ROOTSYS=/usr/local
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib/root:
    if [ -z ${ROOTSYS} ] ; then
	echo "Error: ROOT environment cannot be configured!"
    else echo "Setup root enviroment for user."
    fi
fi

if [[ `root-config --platform` == "macosx" ]]; then

    # With Fink ROOT installations, DYLD_LIBRARY_PATH doesn't have
    # to be defined for ROOT to work. So let's leave the test for it...
    export DYLD_LIBRARY_PATH=${ANALYZER_LIB_PATH}:${DYLD_LIBRARY_PATH}
    
else    
    
    if [ ! $LD_LIBRARY_PATH ]; then
        echo "Warning: so far you haven't setup your ROOT enviroment properly (no LD_LIBRARY_PATH): FrameWork will not work"
    fi
    export LD_TMP_LIBRARY_PATH=${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=${ANALYZER_LIB_PATH}:${LD_LIBRARY_PATH}

fi


export PATH=${ANALYZER_BIN_PATH}:${PATH}
export PYTHONPATH=${ANALYZER_DIR}/python:${PYTHONPATH}
export PAR_PATH=./:${ANALYZER_LIB_PATH}

python ${ANALYZER_DIR}/python/local_check.py

if [ ! -d ${ANALYZER_LOG_PATH} ]; then
    echo Directory ${ANALYZER_LOG_PATH} does not exist ... creating it
    mkdir ${ANALYZER_LOG_PATH}
fi

echo "Running analysis from" $HOSTNAME " in directory: " $PWD

#clean up all emacs tmp files
#clean_emacs