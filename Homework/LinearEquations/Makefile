CFLAGS 		  +=  $(shell gsl-config --cflags)
LDLIBS 		  +=  $(shell gsl-config --libs)
LDLIBS 		  +=  -lm -Wall

CC  		 =  gcc
PIPE     	 =  cat

EXE   	   	 =  main
OBJ  		 =  main.o gramSchmidt.o backsub.o extraFuncs.o timer.o runtime.o
HDRS    	 =  gramSchmidt.h backsub.h extraFuncs.h timer.h
TXT     	 =  out.txt
PLT              =  partC.png
PLTGSL           =  partC_gsl.png
PLTDATA          =  QRtimer.txt
PLTDATAGSL       =  QRtimer_gsl.txt

.PHONEY: default
default: print $(PLT) $(PLTGSL)

print: $(TXT)
	$(PIPE) $<

partC.png: partC.gpi $(PLTDATA)
	pyxplot $<

partC_gsl.png: partC_gsl.gpi $(PLTDATAGSL)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(PLTDATA) $(PLTDATAGSL) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) $(PLTDATAGSL) $(PLTDATA) $(PLT) $(PLTGSL)

