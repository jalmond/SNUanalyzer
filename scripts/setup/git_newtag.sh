if [[ $ANALYZER_DIR == "" ]];
    then
    source $ANALYZER_DIR/setup.sh
fi

########## Tag index
itag=".1"
tagname=$CATVERSION$itag


diff  $ANALYZER_DIR/setup.sh $ANALYZER_DIR/scripts/setup/tag_setup.sh >> SetupCheck.txt

setup_is_different="False"
while read line
do
    setup_is_different="True"
done  < SetupCheck.txt

rm SetupCheck.txt
if [[ $setup_is_different == "True" ]];
then
    echo "setup.sh is changed. Make changes to scripts/setup/tag_setup.sh"
    return
fi

notnew_tag=False
while read line
do
    if [[ $line  == *$tagname* ]]; then
	notnew_tag=True
    fi	
done  < /data1/SNUAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt

if [[ $notnew_tag == "False" ]];then
    echo "$tagname (HEAD)" >> LatestTag80X.txt
    
    while read line
    do
	if [[ $line != *"HEAD"* ]];then
	    echo "$line" >> LatestTag80X.txt
	else
	    suffix="(HEAD)"
	    sline=${line%$suffix}
	    echo "$sline" >> LatestTag80X.txt
	fi
    done  < /data1/SNUAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt
    mv LatestTag80X.txt /data1/SNUAnalyzer_rootfiles_for_analysis/CATTag/LatestTag80X.txt
else
    echo "Not adding a new tag"
fi

rm $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "export CATVERSION="$CATVERSION >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "### If there is a small bug/new code then new subtag is made"  >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "export tag_numerator='"$itag"'"  >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "if [[ \$1 == '"branch"' ]];"  >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    then" >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    export SNUTAG=" >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "else" >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "    export SNUTAG=$CATVERSION"$itag >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "fi" >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh
echo "" >> $ANALYZER_DIR/scripts/setup/SetBrachAndTag.sh

sendemail=false

tag_exists=false
git tag>> gitcheck.txt
while read line
  do
  if [[ $line  == *$tagname* ]];
      then
      tag_exists=true
  fi
done < gitcheck.txt

echo $tag_exists
rm gitcheck.txt


cp $ANALYZER_DIR/scripts/setup/tag_setup.sh $ANALYZER_DIR/setup.sh
git commit -a

#"New Tag: "$tagname  

if [[ $1 == "" ]];
    then
    echo "Making new tag "$tagname
fi



if [[ $tag_exists == "true" ]];
    then
    echo "Deleting branch "$tagname
    git tag -d $tagname
    git push origin :refs/tags/$tagname
    
fi

if [[ $1 == "-e" ]];
    then
    sendemail=true;
fi


git tag $tagname
git push --tags

if [[ $sendemail == "true" ]]; then
    python $ANALYZER_DIR/python/NewTagEmail.py -t $tagname
fi