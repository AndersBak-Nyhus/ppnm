CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

OUTPUT  =  box_solutions.txt my_timing.txt gsl_timing.txt
EXE   	=  main
OBJ  	=  main.o jacobi.o utilities.o diffClock.o test_runtime.o
HDRS    =  jacobi.h utilities.h diffClock.h
IMGFLS   =  plot.png part_C_gsl.png part_C_myMethod.png
TXT     =  output.txt
PLT          =  part_C_myMethod.png
PLTGSL       =  part_C_gsl.png
PLTDATA      =  QRclock_my.txt
PLTDATAGSL   =  QRclock_gsl.txt

.PHONEY: default
default: print $(IMGFLS)

print: $(TXT)
	$(PIPE) $<

plot.png: plot.gpi $(OUTPUT)
	pyxplot $<

part_C_myMethod.png: part_C_myMethod.gpi my_timing.txt
	pyxplot $<

part_C_gsl.png: part_C_gsl.gpi gsl_timing.txt
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)
