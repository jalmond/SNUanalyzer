 # $Id: Makefile 120 2013-11-26 16:42:57Z jalmond $
 ###########################################################################
 # @Project: SFrame - ROOT-based analysis framework for ATLAS              #
 # @Package: User                                                          #
 #                                                                         #
 # @author John       <John.Almond@cern.ch>                - SNU           #
 #                                                                         #
 # Makefile compiling the SNUTree    library.                              #
 #                                                                         #
 ###########################################################################
# Package information
LIBRARY = SKTree
OBJDIR  = $(OBJ)
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include


# Override default target so we can install the rootmap file
SKTreedefault: default rootmap SKTreeUnit

# Include the generic compilation rules
include $(ANALYZER_DIR)/Makefile.common

SKTreeUnit: SKTreeUnit.o $(SHLIBFILE)
	@echo "Linking " $@
	@$(LD) $(LDFLAGS)  -O2 $(OBJDIR)/SKTreeUnit.o -L$(ANALYZER_LIB_PATH) \
                $(ROOTLIBS) -lSKTree -o $@

SKTreeUnit.o: SKTreeUnit.cc include/KParticle.h
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -O2 -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)

#install the rootmap file
rootmap: $(ANALYZER_DIR)/${Flag}Core/SKTree/util/SKTree.rootmap
	@echo "Installing rootmap " $<
	(cp $< $(ANALYZER_LIB_PATH)/.)



