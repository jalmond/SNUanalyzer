snuagpath=/data1/${Flag}Analyzer_rootfiles_for_analysis/SNUTag/
if [ $HOSTNAME == "cmscluster.snu.ac.kr" ];
then
    tagpath=/data4/LocalNtuples/${Flag}Analyzer_rootfiles_for_analysis/SNUTag/
else
    tagpath=/data1/${Flag}Analyzer_rootfiles_for_analysis/SNUTag/
fi

if [[ $USER == "jalmond" ]];
    then
    if [[ $SNUTAG == "" ]];
	then
	echo "This is a new git checkout of the main branch. You need the tag setup to run code" 
	cp $ANALYZER_DIR/scripts/setup/tag_setup.sh $ANALYZER_DIR/setup.sh
	export SNUTAG=$SNUVERSION$tag_numerator
    fi
fi
if [[ $SNUTAG == "" ]];
    then
    echo "You are not running from a tag. You are running from Branch SnuAnalyzer_13TeV"
    exit 1
fi

if [[ $SNUTAG == "v7-6-3.2" ]];
    then
    echo "You are running on a tag with a bug in pileup weighting. update tag"
    exit 1
fi

latest_tag=""
while read line
do
    if [[ $line == *"HEAD"* ]];
    then
	sline=$(echo $line | head -n1 | awk '{print $1}')
	latest_tag=$sline
    fi
done < $tagpath/LatestTag94X.txt


if [[ $latest_tag == $SNUTAG ]];then
    
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    echo "Using latest tag "$SNUTAG
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
else
    
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    echo "Newer SNUAnalzer tag available: "
    echo "Current tag "$SNUTAG
    echo "Latest tag is "$latest_tag
    echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

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
    done < $tagpath/LatestTag94X.txt


    for ntag in  ${NEWTAGS[@]};
    do
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	echo "Tag: " $ntag  "(summary of changes wrt previous tag)"
	echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	
	while read line
	do
	    echo $line
	done < $tagpath/TagDiff_${ntag}.txt
    done
fi