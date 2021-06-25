CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o minimization.o
HD      =  minimization.h
TXT     =  output.txt
DATA    =  data.txt

.PHONEY: default
default: print


print: $(TXT)

$(TXT) : $(EXE)
	./$< $(DATA) 1> $(TXT)

$(EXE): $(OBJ) $(HD)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)
