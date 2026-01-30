rootlibs:=$(shell root-config --libs)
rootflags:=$(shell root-config --cflags)
rootflags2 := -I${ROOTSYS}/include
corlibincl := $(shell correction config --incdir)
#corliblib := /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/lib
corliblib = $(shell a=`correction config --libdir`; echo $${a%%python*lib})
corlibflag := $(shell correction config --cflags --ldflags --rpath)
SOFLAGS = -shared

LD = g++ -m64 -g -Wall

SRCDIR=src
COMMONDIR=../CommonTools/src

CXXFLAGS = -O0 -g -Wall -fmessage-length=0 $(rootflags) -fpermissive -fPIC -pthread -DSTANDALONE -I. -I$(corlibincl) -I$(COMMONDIR) -I$(SRCDIR)

# source files
SRCS_SRC     := $(wildcard $(SRCDIR)/*.cpp)
SRCS_COMMON  := $(wildcard $(COMMONDIR)/*.cpp)
#SRCS_COMMON  := $(filter-out $(SRC_COMMON)/rootdict.C \
#	                         $(SRC_COMMON)/JetMETObjects_dict.C, \
#							 $(wildcard $(SRC_COMMON)/*.cpp))

# object files (same path, .cpp -> .o)
OBJS_SRC     := $(patsubst %.cpp,%.o,$(SRCS_SRC))
OBJS_COMMON  := $(patsubst %.cpp,%.o,$(SRCS_COMMON))

OBJS := $(OBJS_SRC) $(OBJS_COMMON) \
        $(COMMONDIR)/JetMETObjects_dict.o \
        $(COMMONDIR)/rootdict.o

LIBS_EXE = $(rootlibs) -lMathMore -lGenVector -L$(corliblib) -lcorrectionlib
LIBS = $(rootlibs)

TARGET =	nanoaodrdataframe

all:	$(TARGET) libnanoadrdframe.so 

clean:
	rm -f $(OBJS) $(TARGET) libnanoaodrdframe.so $(COMMONDIR)/JetMETObjects_dict.C $(COMMONDIR)/rootdict.C JetMETObjects_dict_rdict.pcm rootdict_rdict.pcm

$(COMMONDIR)/rootdict.C: $(COMMONDIR)/NanoAODAnalyzerrdframe.h $(SRCDIR)/TopLFVAnalyzer.h $(COMMONDIR)/TauFakeFactorAnalyzer.h $(SRCDIR)/SkimEvents.h $(COMMONDIR)/Linkdef.h
	rm -f $@
	rootcint $@ -I$(corlibincl) -I$(COMMONDIR) -I$(SRCDIR) $^
	rm -f rootdict_rdict.pcm
	ln -s $(COMMONDIR)/rootdict_rdict.pcm .

	
libnanoadrdframe.so: $(OBJS)
	$(LD) $(SOFLAGS) $(LIBS_EXE) -D_GLIBCXX_USE_CXX11_ABI=0 -o $@ $^


$(COMMONDIR)/JetMETObjects_dict.C: $(COMMONDIR)/JetCorrectorParameters.h $(COMMONDIR)/SimpleJetCorrector.h $(COMMONDIR)/FactorizedJetCorrector.h $(COMMONDIR)/JetResolutionObject.h $(COMMONDIR)/LinkdefJetmet.h
	rm -f $@
	$(ROOTSYS)/bin/rootcint -f $@ -c $(rootflags2) $^
	ln -s $(COMMONDIR)/JetMETObjects_dict_rdict.pcm .

$(COMMONDIR)/rootdict.o: $(COMMONDIR)/rootdict.C
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(COMMONDIR)/JetMETObjects_dict.o: $(COMMONDIR)/JetMETObjects_dict.C
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(COMMONDIR)/%.o: $(COMMONDIR)/%.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS_EXE)
