# $Id: Makefile
###########################################################################                                                                                                    
# @Project: LQAnalyzer - ROOT-based analysis framework for Korea CMS      #                                                                                                    
#                                                                         #                                                                                                    
# @author John Almond       jalmond@cern.ch>           - SNU              #                                                                                                    
# Top level Makefile for compiling all the LQAnalyzer code                #                                                                                                    
#                                                                         #                                                                                                    
###########################################################################                                                                                                    


# Package information
LIBRARY = Ntuplecore
OBJDIR  = $(OBJ)
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

coredefault: default $(ANALYZER_BIN_PATH)/run_lqanalyzer   

# Include the generic compilation rules
include $(ANALYZER_DIR)/Makefile.common

INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/SKTree/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Cycle/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Plotting/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Selection/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/AnalysisCore/include/


$(ANALYZER_BIN_PATH)/run_lqanalyzer: run_lqanalyzer.o $(SHLIBFILE)                                                                           
	@echo "Linking " $(LD) $(LDFLAGS) $@                                                                                                                     
	@$(LD) $(LDFLAGS) $(OBJDIR)/run_lqanalyzer.o -L$(ANALYZER_LIB_PATH) -lNtuplecore -lSKTree -lHist -lAnalysisCore   \
		$(ROOTLIBS) -lPyROOT  -o $@ 

run_lqanalyzer.o: app/run_lqanalyzer.cxx                                                                                                       
	@echo "Compiling $<"                                                                                                                    
	@mkdir -p $(OBJDIR)                                                                                                                     
	@$(CXX) $(CXXFLAGS) -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)