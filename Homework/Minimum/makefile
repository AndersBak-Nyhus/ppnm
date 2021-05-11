CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o minimization.o
HDRS    =  minimization.h
TXT     =  output.txt
IMGFL   =  plot.png
GPI     =  plot.gpi
DATA    =  data.txt

.PHONEY: default
default: print #$(IMGFL)

#hydrogen_plot.png: $(HGPI) $(HDATA)
#	pyxplot $<

#convergence_plot.png: $(CGPI) $(CDATA)
#	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(DATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)
