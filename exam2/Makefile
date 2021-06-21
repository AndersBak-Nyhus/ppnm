CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

TXT = output.txt
IN = input.txt
HD = input2Array.h bilinear.h binsearch.h
OBJ = input2Array.o main.o binsearch.o bilinear.o
EXE = main

all: output.txt


$(TXT): main 
	./$<  >$(TXT)

$(EXE): $(OBJ)  $(HD)


.PHONEY: clean
clean:
	$(RM) *.o out* *.svg test* log* main
