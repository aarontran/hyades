# Super ugly
CC= DYLD_LIBRARY_PATH=/Users/atran/opt/miniconda3/envs/wham/lib clang++

INCL= -I/Users/atran/opt/miniconda3/envs/wham/include
LIBS= -L/Users/atran/opt/miniconda3/envs/wham/lib

FLAGS= -Wall -fopenmp -g
#FLAGS= -lgomp -fopenmp
#FLAGS= -fopenmp -Ofast -mcpu=native -mtune=native -fno-strict-aliasing


# Quick reference.
# $@ = file name of the target
# $^ = names of all the prerequisites, with spaces between them

%.o: src/%.cc Makefile
	$(CC) $(INCL) $(FLAGS) -c $<

hyades: field.o field_advance.o field_constructor.o hyades.o interp.o particle.o particle_move.o random.o
	$(CC) $(INCL) $(LIBS) $(FLAGS) -o hyades $^

all: hyades

clean: 
	rm -f hyades
	rm -f *.o
