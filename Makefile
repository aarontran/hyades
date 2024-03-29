#USER=input/pcai.cc
USER=input/column-nz32.cc

# use to select particle shape
# terrible hack need to change
#FLAGS_USER=-DSHAPE_NGP
FLAGS_USER=-DSHAPE_QS
#FLAGS_USER=-DSHAPE_CIC
#FLAGS_USER=-DSHAPE_TSC

# Personal Macbook
#CC= DYLD_LIBRARY_PATH=/Users/atran/opt/miniconda3/envs/wham/lib clang++
CC= DYLD_LIBRARY_PATH=/Users/atran/opt/miniconda3/envs/wham/lib h5c++
INCL= -I/Users/atran/opt/miniconda3/envs/wham/include
LIBS= -L/Users/atran/opt/miniconda3/envs/wham/lib
#INCL=
#LIBS=
#FLAGS= -Wall -fopenmp -g -O0 -std=c++17
#FLAGS= -lgomp -fopenmp
FLAGS= -fopenmp -Ofast -mcpu=native -mtune=native -std=c++17
#-fno-strict-aliasing  # this degrades performance by 20% ??
# need std c++17 to get std::filesystem library

## AOCC (on Purdue Anvil using module aocc/3.1.0)
#CC=mpicxx
#INCL=
#LIBS=
#FLAGS= -Ofast -march=znver3 -fno-strict-aliasing

## GCC
#CC=g++
#INCL=
#LIBS=
#FLAGS= -O2 -fno-strict-aliasing -march=native

## Intel (on Purdue Anvil using modules intel/19.0.5.281 impi/2019.5.281 hdf5/1.10.7)
#CC=h5c++
#INCL=
#LIBS=
#FLAGS= -O2 -march=core-avx2 -ipo
##-qopt-zmm-usage=high 	# for intel architectures
##-xCORE-AVX512


# Quick reference.
# $@ = file name of the target
# $^ = names of all the prerequisites, with spaces between them

hyades: field.o field_advance.o field_constructor.o field_dump.o hyades.o interp.o particle.o particle_move.o particle_dump.o random.o timer.o $(addsuffix .o, $(basename $(notdir $(USER))))
	$(CC) $(INCL) $(LIBS) $(FLAGS) $(FLAGS_USER) -o hyades $^

%.o: src/%.cc Makefile
	$(CC) $(INCL) $(FLAGS) $(FLAGS_USER) -c $<

$(addsuffix .o, $(basename $(notdir $(USER)))): $(USER) Makefile
	$(CC) $(INCL) $(FLAGS) $(FLAGS_USER) -c $<

all: hyades

clean: 
	rm -f   hyades
	rm -f   *.o
	rm -rf  output
