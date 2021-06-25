CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

OUTPUT  =  box.txt timing.txt gsl_timing.txt
EXE   	=  main
OBJ  	=  main.o jacobi.o extraFuncs.o timer.o runtime.o
HDRS    =  jacobi.h extraFuncs.h timer.h
IMG   =  plot.png part_C_gsl.png part_C.png
TXT     =  output.txt
PLT          =  part_C.png
PLTGSL       =  part_C_gsl.png
PLTDATA      =  QRclock_my.txt
PLTDATAGSL   =  QRclock_gsl.txt

.PHONEY: default
default: print $(IMG)

print: $(TXT)
	$(PIPE) $<

plot.png: plot.gpi $(OUTPUT)
	pyxplot $<

part_C.png: part_C.gpi timing.txt
	pyxplot $<

part_C_gsl.png: part_C_gsl.gpi gsl_timing.txt
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMG) *.png
