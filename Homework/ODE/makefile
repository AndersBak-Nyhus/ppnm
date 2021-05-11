CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

OUTPUT  =  harmonic_ODE_output.txt SIR_ODE_output.txt SIR2_ODE_output.txt threeprob_ODE_output.txt
EXE   	=  main
OBJ  	=  main.o RK-ODE.o functions.o
HDRS    =  RK-ODE.h functions.h
IMGFL   =  harmonic_plot.png sir_plot.png sir2_plot.png threeprob_plot.png
TXT     =  output.txt

.PHONEY: default
default: print $(IMGFL)

print: $(TXT)
	$(PIPE) $<

harmonic_plot.png: harmonic_plot.gpi $(HARMOUTPUT)
	pyxplot $<

sir_plot.png: sir_plot.gpi $(SIROUTPUT)
	pyxplot $<

sir2_plot.png: sir2_plot.gpi $(SIR2OUTPUT)
	pyxplot $<

threeprob_plot.png: threeprob_plot.gpi $(TPROBOUTPUT)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
