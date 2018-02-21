##################################################################################
## Author: P T Surukuchi
## Date: June 04, 2016
## OscSens Makefile
## TODO: Force recompilation when header files change
## 
##################################################################################

## Avoids trouble on systems where the SHELL variable is inherited from the environment
SHELL = /bin/sh

## Clear suffix lists
.SUFFIXES:
## Apply implicit rules on the following suffixes
.SUFFIXES: .cc .o

## C++ flags to be used, 
## O3: Extreme optimization (performance and code size) at expense of compile time
##	and debuggability. Can also use O2 or O1 based on what you want
## --std=c++11 : Use c++11 
## -fPIC: Generate position-independent code (PIC) suitable for use in a shared library
## `root-config --cflags` : Flags for compiling a source file against the ROOT header files.
## -Wpedantic : Issue all the warnings demanded by strict ISO C and ISO C++.
## -Wno-vla-extension : Suppress warnings related to variable length arrays (VLAs are supported for C99 and over)
## -Wall, -Wextra : Enforces strict compile-time warnings
## -g : Request gcc to emit extra information for debugging
##	 In case of debugging, replace -O3 with -Og for optimizing debugging
CXXFLAGS = -O3 --std=c++11 -fPIC `root-config --cflags` -Wpedantic -Wno-vla-extension -Wall -Wextra

## root-config --libs : Use root-related libraries
LDFLAGS = `root-config --libs` -L./

## Include locations to search for source files
INC = -I./ -I./src -I./src/DataExtractionCode

## Search in the source directory for .cc and .hh files
VPATH = ./src

## Search in the data extraction subdirectory for .cc and .hh files
VDATAEXTPATH =  ./src/DataExtractionCode

## Directory containing .o files
OBJDIR = ./obj

## Directory containing the library
LIBDIR = ./lib

## Directory containing code for running some basic tests
TESTSDIR = ./tests/

## Define names of the objects (saved in ./obj) of files used by the executables
OBJS = $(addprefix $(OBJDIR)/, TDetector.o TDetectorLocation.o TReactor.o TExperiment.o TOscillationSimulator.o TOscillationModelBuilder.o TMinimization.o TMacroExtractor.o TCovarianceMatrixInterface.o TMacroInterface.o TMCToyInterface.o TThrowMCToy.o TDataExtractor.o TCovarianceMatrixGenerator.o TCovUtilities.o TDetectorResponse.o TDetectorResponseInterface.o TOscSensUtils.o TFileUtilities.o)

## Define oscillation Sensitivity library
LIBOSC = $(LIBDIR)/libOscSens.so

## Generate library using the object files
$(LIBOSC): $(OBJS)
	ar rs $@ $(OBJS)

## Create .o files in the object directory
$(OBJDIR)/%.o : $(VPATH)/%.cc $(VPATH)/%.hh
	$(CXX) $(CXXFLAGS) $< -c -o $@

## Create .o files in the object directory
$(OBJDIR)/%.o : $(VDATAEXTPATH)/%.cc $(VDATAEXTPATH)/%.hh
	$(CXX) $(CXXFLAGS) $(INC) $< -c -o $@

## All the executable files generated separately ($ make executable) 
% : %.cc $(LIBOSC)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LDFLAGS) $(LIBOSC) -o $@

## All the executable files generated separately ($ make executable) 
% : ./aux/%.cc $(LIBOSC)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LDFLAGS) $(LIBOSC) -o $@

## All the executable files generated separately ($ make executable) 
% : ./aux/ToyCovExecs/%.cc $(LIBOSC)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LDFLAGS) $(LIBOSC) -o $@

## All the executable files generated separately ($ make executable) 
% : $(TESTSDIR)%.cc $(LIBOSC)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LDFLAGS) $(LIBOSC) -o $@

## Make sure to check if obj directory exists and create if it doesn’t
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)

## Make sure to check if lib directory exists and create if it doesn’t
$(LIBOSC): | $(LIBDIR)
$(LIBDIR):
	mkdir $(LIBDIR)

## Remove .o, .so files
.PHONY: clean
clean:
	-rm -rf $(OBJS) $(LIBOSC) $(TESTSDIR)/*.mac $(TESTSDIR)/*.root
