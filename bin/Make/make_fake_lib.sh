cd $ANALYZER_DIR/SNUAnalysis/AnalyzerTools/HNCommonLeptonFakes/conf/
make -f Makefile.StandAlone
cd -
cd $ANALYZER_LIB_PATH
if [[ -f libHNCommonLeptonFakes.so ]];
    then
    rm libHNCommonLeptonFakes.so
fi
cp $ANALYZER_DIR/SNUAnalysis/AnalyzerTools/HNCommonLeptonFakes/Root/libHNCommonLeptonFakes.so .
cd $ANALYZER_DIR/SNUAnalysis/AnalyzerTools/
make distclean
make
cd $ANALYZER_DIR