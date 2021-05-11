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
IMGFL   =  hydrogen_plot.png convergence_plot.png
HGPI    =  hydrogen_plot.gpi
CGPI    =  convergence_plot.gpi
HDATA   =  hydrogen.txt
CDATA   =  convergence.txt

.PHONEY: default
default: print $(IMGFL)

hydrogen_plot.png: $(HGPI) $(HDATA)
	pyxplot $<

convergence_plot.png: $(CGPI) $(CDATA)
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(HDATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)
