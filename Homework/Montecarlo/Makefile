CC = gcc
CFLAGS =-std=gnu99 -O -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

TXT = out.txt
IMG = error.png
EXE = main
OBJ = Montecarlo.o main.o

default: out.txt error.png
	cat out.txt

error.png: error.pgi error.txt
	pyxplot $<

$(TXT): $(EXE)
	./$(EXE) > $(TXT)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE): $(OBJ)
main.c: Montecarlo.c
clean:
	$(RM) $(OBJ) $(TXT) $(EXE) $(IMG) *.txt
