ObjSuf = o
SrcSuf = cc
ExeSuf =
DllSuf = so
OutPutOpt = -o
HeadSuf = h

ROOTCFLAGS = $(shell root-config --cflags) -O2
ROOTLIBS = $(shell root-config --libs) -O2

# Linux with egcs

CXX = g++ -std=c++11
CXXFLAGS += $(ROOTCFLAGS) -I./src/ -Wall 

LD = g++ -std=c++11
LDFLAGS_LOCAL += $(ROOTLIBS) -Wall -l TreePlayer -lboost_program_options -Wno-parentheses -Wno-unused-but-set-variable -L/nfs_scratch/jmadhusu/CMSSW_10_2_10/lib/slc7_amd64_gcc700 -lHTTutilitiesJet2TauFakes
LDFLAGS_CONDOR+= $(ROOTLIBS) -Wall -l TreePlayer -lboost_program_options -Wno-parentheses -Wno-unused-but-set-variable -L$$PWD -lHTTutilitiesJet2TauFakes
LDPATHCONDOR = -L$$PWD

SOFLAGS = -shared
LIBS =

SRCDIR = src
OBJDIR = obj

#------------------------------------------------------------------------------
SOURCES = $(wildcard $(SRCDIR)/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
#-----------------------------------------------------------------------------

all: tautau_analyzer

Analyzer_tautau: $(OBJECTS)  obj/tautau_analyzer.O
	@$(LD) -o $@ $^ $(LDFLAGS_LOCAL) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@$(CXX) -o $@  -c $< $(CXXFLAGS)

$(OBJDIR)/%.O: $(SRCDIR)/%.C
	@$(CXX) -o $@  -c $< $(CXXFLAGS)

condor: $(OBJECTS)  obj/Analyzer_tautau.O
	@$(LD) -o $@ $^ $(LDFLAGS_CONDOR) $(LIBS)

clean:
	@echo "Cleaning..."
	@ls $(OBJDIR)
	@rm -f $(OBJECTS) obj/*.O

.SUFFIXES: .$(SrcSuf) .cc .o .so
