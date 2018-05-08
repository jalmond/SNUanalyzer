# Package information
LIBRARY = Validation
OBJDIR  = $(OBJ)
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# Overwrite the default rule defined in Makefile.common                                                                                                                      
# Include the generic compilation rules
include $(ANALYZER_DIR)/Makefile.common

INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/SKTree/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Ntuplecore/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Selection/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Plotting/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/AnalysisCore/include/
INCLUDES += -I${BTAGDIR}
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Analysis/AnalyzerTools/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Analysis/AnalyzerTools/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Analysis/Analyzers/include/



ifeq ($(isSLC),False)
	CXXFLAGS= -O2 -std=c++0x  -Wall -fPIC
	CXXFLAGS     += $(ROOTCFLAGS)
endif