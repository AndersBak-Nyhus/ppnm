CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o rootfinding.o gramSchmidt.o backsub.o RK-ODE.o
HDRS    =  rootfinding.h gramSchmidt.h backsub.h RK-ODE.h
TXT     =  output.txt
IMG	=  hydrogen.png convergence.png
HGPI    =  hydrogen_plot.gpi
CGPI    =  convergence_plot.gpi
HDATA   =  hydrogen.txt
CDATA   =  convergence.txt

.PHONEY: default
default: print $(IMG)

hydrogen.png: $(HGPI) $(HDATA)
	pyxplot $<

convergence.png: $(CGPI) $(CDATA)
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(HDATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) $(IMG)
