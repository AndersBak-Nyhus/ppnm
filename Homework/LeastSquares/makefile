CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

INPUT   =  data.txt
OUTPUT  =  fit.plt
EXE   	=  main
OBJ  	=  main.o gramSchmidt.o backsub.o input2Array.o leastSquares.o utilities.o
HDRS    =  gramSchmidt.h backsub.h input2Array.h leastSquares.h utilities.h
IMGFL   =  fit.png
TXT     =  output.txt

.PHONEY: default
default: print $(IMGFL)

print: $(TXT)
	$(PIPE) $<

$(IMGFL): fit.gpi $(INPUT) $(OUTPUT)
	pyxplot $<

$(TXT) : $(EXE) $(INPUT)
	./$< $(INPUT) $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
