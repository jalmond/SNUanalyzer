cd $ANALYZER_DIR/${Flag}Analysis/AnalyzerTools/HNCommonLeptonFakes/conf/; 
make clean -f Makefile.StandAlone; 
if [[ -f ${ANALYZER_LIB_PATH}/libHNCommonLeptonFakes.so ]];
    then
    rm ${ANALYZER_LIB_PATH}/libHNCommonLeptonFakes.so
fi
cd $ANALYZER_DIR