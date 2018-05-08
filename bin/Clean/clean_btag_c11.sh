cd $ANALYZER_DIR/${Flag}Analysis/AnalyzerTools/BTag/BTagC11/conf/; 
make clean -f Makefile.StandAlone; 
if [[ -f ${ANALYZER_LIB_PATH}/libBTagSFUtil.so  ]];
    then
    rm ${ANALYZER_LIB_PATH}/libBTagSFUtil.so
fi
cd $ANALYZER_DIR