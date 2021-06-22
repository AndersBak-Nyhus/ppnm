CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

ERRPLOT =  error_scaling.png
ERRDATA =  error_scaling.txt
ERRGPI  =  error_scaling.gpi

EXE   	=  main
OBJ  	=  main.o montecarlo.o
HDRS    =  montecarlo.h
TXT     =  output.txt

.PHONEY: default
default: print $(ERRPLOT)

$(ERRPLOT): $(ERRGPI) $(ERRDATA)
	pyxplot $<

print: $(TXT)
	$(PIPE) $<

$(TXT) : $(EXE)
	./$< 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)
