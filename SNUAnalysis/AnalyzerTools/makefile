# Package information
LIBRARY = AnalyzerTools
OBJDIR  = $(OBJ)
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# Overwrite the default rule defined in Makefile.common                                                                                                                      
# Include the generic compilation rules
include $(ANALYZER_DIR)/Makefile.common
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Ntuplecore/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/Selection/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/SKTree/include/
INCLUDES += -I$(ANALYZER_DIR)/${Flag}Core/AnalysisCore/include/
INCLUDES += -I${ROCHDIR}
#INCLUDES += -I$(LQANALYZER_DIR)

ifeq ($(isSLC),False)
	CXXFLAGS= -O2 -std=c++0x  -Wall -fPIC
	CXXFLAGS     += $(ROOTCFLAGS)
endif