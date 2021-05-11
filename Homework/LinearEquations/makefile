CFLAGS 		  +=  $(shell gsl-config --cflags)
LDLIBS 		  +=  $(shell gsl-config --libs)
LDLIBS 		  +=  -lm -Wall

CC		  		 =  gcc
PIPE     		 =  cat

EXE   	   	 =  main
OBJ  	  	 	 =  main.o gramSchmidt.o backsub.o utilities.o diffClock.o test_runtime.o
HDRS    	   =  gramSchmidt.h backsub.h utilities.h diffClock.h
TXT     	 	 =  output.txt
PLT          =  part_C_myMethod.png
PLTGSL       =  part_C_gsl.png
PLTDATA      =  QRclock_my.txt
PLTDATAGSL   =  QRclock_gsl.txt

.PHONEY: default
default: print $(PLT) $(PLTGSL)

print: $(TXT)
	$(PIPE) $<

part_C_myMethod.png: part_C_myMethod.gpi $(PLTDATA)
	pyxplot $<

part_C_gsl.png: part_C_gsl.gpi $(PLTDATAGSL)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(PLTDATA) $(PLTDATAGSL) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) $(PLTDATAGSL) $(PLTDATA) $(PLT) $(PLTGSL)
