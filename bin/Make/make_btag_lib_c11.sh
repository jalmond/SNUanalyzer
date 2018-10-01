cd $BTAGDIR/conf
make -f Makefile.StandAlone
cd -
cd $ANALYZER_LIB_PATH
if [[ -f libBTagSFUtil.so ]] ;
    then
    rm libBTagSFUtil.so
fi
cp $BTAGDIR/Root/libBTagSFUtil.so libBTagSFUtil.so

cd $ANALYZER_DIR/SNUAnalysis/AnalyzerTools/

if [[ $1 == "False" ]]; then
    make distclean
    make
fi
cd $ANALYZER_DIR/

