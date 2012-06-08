# To compile your files simply type `make'
# NOTE: files must have the .C or .cpp suffix
# Filenames must be included in the EXS variable

include ../../makevars

# Set these paths to the location of the GA library and headers.
#GA_INC_DIR= /usr/local/include
#GA_LIB_DIR= /usr/local/lib
GA_INC_DIR= ../..
GA_LIB_DIR= ../../ga

INC_DIRS= -I$(GA_INC_DIR)
LIB_DIRS= -L$(GA_LIB_DIR)

# if you need to compile more files simply insert their name 
# in the EXS variable (NOTE: do not use suffices), in example
# to compile the files ex1.cpp (or .C), foo1.cpp, etc. use: 
# EXS=ex1 foo1 ex3\
# ex4 foo2
EXS=nfiles_ndisks

.SUFFIXES: .C .cpp
.C.o .cpp.o:
	$(CXX) $(CXXFLAGS) $(INC_DIRS) -c $<

all: $(EXS)

# Use this for non-gnu make
#$(EXS): $$@.o
#	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm $(CXX_LIBS)

# Use this for gnu make
$(EXS): %: %.o
	$(CXX) $@.o -o $@ $(LIB_DIRS) -lga -lm $(CXX_LIBS)

clean:
	$(RM) $(EXS)
	$(RM) *.o *~ *.bak *.pixie core
	$(RM) test_results.txt test_stats.txt
	$(RM) *.dat
	$(RM) *.out *.exe vc* *.pdb
