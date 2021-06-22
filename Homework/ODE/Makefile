CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

OUTPUT  =  harmonic_ODE_output.txt SIR_ODE_output.txt SIR2_ODE_output.txt threeprob_ODE_output.txt
EXE   	=  main
OBJ  	=  main.o RK-ODE.o functions.o
HDRS    =  RK-ODE.h functions.h
IMG   =  harmonic.png SIR.png SIR2.png threebody.png
TXT     =  output.txt

.PHONEY: default
default: print $(IMG)

print: $(TXT)
	$(PIPE) $<

harmonic.png: harmonic_plot.gpi $(HARMOUTPUT)
	pyxplot $<

SIR.png: sir_plot.gpi $(SIROUTPUT)
	pyxplot $<

SIR2.png: sir2_plot.gpi $(SIR2OUTPUT)
	pyxplot $<

threebody.png: threeprob_plot.gpi $(TPROBOUTPUT)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMG)
