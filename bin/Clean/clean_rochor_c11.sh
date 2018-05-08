cd $ANALYZER_DIR/${Flag}Analysis/AnalyzerTools/rochcor2016/conf/;
make clean -f Makefile.StandAlone; 
if [[ -f ${ANALYZER_LIB_PATH}/librochcor2016.so  ]];
    then
    rm ${ANALYZER_LIB_PATH}/librochcor2016.so
fi
cd $ANALYZER_DIR