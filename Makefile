 # $Id: Makefile
 ###########################################################################
 # @Project: LQAnalyzer - ROOT-based analysis framework for Korea CMS      #
 #                                                                         #
 # @author John Almond       jalmond@cern.ch>           - SNU              #
 # Top level Makefile for compiling all the LQAnalyzer code                #
 #                                                                         #
 ###########################################################################

all: sktree core AnalysisCore plotting selection cycle 

sktree::
	(cd SKTree; make)

core::
	(cd LQCore/core; make)

AnalysisCore::
	(cd LQCore/AnalysisCore; make)

plotting::
	(cd LQCore/Plotting; make)

selection::
	(cd LQCore/Selection; make)

cycle::
	(cd LQCycle; make)

clean::
	(cd SKTree; make clean)	
	(cd LQCore/core; make clean)
	(cd LQCore/AnalysisCore; make clean)
	(cd LQCore/Plotting; make clean)
	(cd LQCore/Selection; make clean)
	(cd LQCycle; make clean)

distclean::
	(cd SKTree; make distclean)
	(cd LQCore/core; make distclean)
	(cd LQCore/AnalysisCore; make distclean)
	(cd LQCore/Plotting; make distclean)
	(cd LQCore/Selection; make distclean)
	(cd LQCycle; make distclean)