 # $Id: Makefile
 ###########################################################################
 # @Project: SNUAnalyzer - ROOT-based analysis framework for Korea CMS      #
 #                                                                         #
 # @author John Almond       jalmond@cern.ch>           - SNU              #
 # Top level Makefile for compiling all the SNUAnalyzer code                #
 #                                                                         #
 ###########################################################################

all: tagcheck btag roch fakes sktree AnalysisCore Ntuplecore plotting selection analysis 

sktree::
	(cd SNUCore/SKTree; make)

Ntuplecore::
	(cd SNUCore/Ntuplecore; make)

roch::
	(bash bin/Make/make_rocher_lib.sh; cd ${ANALYZER_DIR} )

AnalysisCore::
	(cd SNUCore/AnalysisCore; make)

plotting::
	(cd SNUCore/Plotting; make)

selection::
	(cd SNUCore/Selection; make)

analysis::
	(cd SNUAnalysis/AnalyzerTools; make)
	(cd SNUAnalysis/Analyzers; make)
	(cd SNUAnalysis/SKTreeMaker; make)
	(cd SNUAnalysis/Validation; make)


fakes::
	(cd ${ANALYZER_DIR}/SNUAnalysis/AnalyzerTools/HNCommonLeptonFakes/conf/; make -f Makefile.StandAlone; cd ${ANALYZER_LIB_PATH} ;rm libHNCommonLeptonFakes.so ; cp ${ANALYZER_DIR}/SNUAnalysis/AnalyzerTools/HNCommonLeptonFakes/Root/libHNCommonLeptonFakes.so .; cd ${ANALYZER_DIR} )

btag::  
	(bash bin/Make/make_btag_lib.sh; cd ${ANALYZER_DIR} )	

tagcheck::
	(source bin/CheckNewTagCompiler.sh ${CHECKTAGFILE})

clean::
	(cd SNUCore/SKTree; make clean)
	(cd SNUCore/Ntuplecore; make clean)
	(cd SNUCore/AnalysisCore; make clean)
	(cd SNUCore/Plotting; make clean)
	(cd SNUCore/Selection; make clean)
	(cd SNUAnalysis/AnalyzerTools; make clean)
	(cd SNUAnalysis/Analyzers; make clean)
	(cd SNUAnalysis/SKTreeMaker; make clean)
	(cd SNUAnalysis/Validation; make clean)
	(bash bin/Clean/clean_fake.sh)
	(bash bin/Clean/clean_rochor.sh)
	(bash bin/Clean/clean_btag.sh)

distclean::
	(cd SNUCore/SKTree; make distclean)
	(cd SNUCore/Ntuplecore; make distclean)
	(cd SNUCore/AnalysisCore; make distclean)
	(cd SNUCore/Plotting; make distclean)
	(cd SNUCore/Selection; make distclean)
	(cd SNUAnalysis/AnalyzerTools; make distclean)
	(cd SNUAnalysis/Analyzers; make distclean)
	(cd SNUAnalysis/SKTreeMaker; make distclean)
	(cd SNUAnalysis/Validation; make distclean)


	(bash bin/Clean/clean_fake.sh)
	(bash bin/Clean/clean_rochor.sh)
	(bash bin/Clean/clean_btag.sh)

