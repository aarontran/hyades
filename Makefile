# Super ugly
CC= DYLD_LIBRARY_PATH=/Users/atran/opt/miniconda3/envs/wham/lib clang++

INCL= -I/Users/atran/opt/miniconda3/envs/wham/include
LIBS= -L/Users/atran/opt/miniconda3/envs/wham/lib
FLAGS= -Wall -fopenmp -g
#FLAGS= -lgomp -fopenmp
#FLAGS= -Ofast -mcpu=native -mtune=native -fno-strict-aliasing


# Quick reference.
# $@ = file name of the target
# $^ = names of all the prerequisites, with spaces between them

%.o: src/%.cc Makefile
	$(CC) $(INCL) $(LIBS) $(FLAGS) -c $<

hyades: fields.o hyades.o particles.o random.o
	$(CC) $(INCL) $(LIBS) $(FLAGS) -o hyades $^

all: hyades

clean: 
	rm hyades
	rm *.o
