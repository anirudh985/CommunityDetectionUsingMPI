# This is a very simple makefile for building the Lisp interpreter
# project when using C++ on stdlinux. Feel free to add any improvements:
# e.g. pattern rules, automatic tracking of dependencies, etc. There
# is a lot of info about "make" on the web.

# C++ compiler
CXX = mpic++

# C++ compiler flags
CXXFLAGS = -Ofast -fopenmp -std=c++0x

# Creating a .o file
COMPILE = $(CXX) $(CXXFLAGS) -c

# Name of the executable; should match the contents of Runfile
EXE = louvainMain

# All object files
OBJS = Graph.o louvainMultiPhase.o louvainMain.o utilityFunctions.o parallelLouvain.o ReadFile.o ReformatFile.o

# The first target is the one that is executed when you invoke
# "make". The line describing the action starts with <TAB>. Variable
# "$@" stands for the current target. 

$(EXE) : $(OBJS) # depends on all object files
	$(CXX) $(CXXFLAGS) $^ -o $@

# An object file is dependent on the corresponding source file

%.o : %.cpp
	$(COMPILE) -o $@ $<
