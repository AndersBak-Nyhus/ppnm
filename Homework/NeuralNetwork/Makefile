CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


TXT = out.txt
IMG = plot.png
PLTXT = Pts.txt Funcs.txt
EXE = main
OBJ = main.o ann.o minimization.o

.PHONEY: default

default: $(TXT) $(IMG)
	cat $(TXT)

$(IMG): plots.pgi $(PLTXT)
	pyxplot $<


$(TXT): $(EXE)
	./$(EXE) > $(TXT)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

main: $(OBJ)
main.c: ann.c
ann.c: minimization.c

.PHONEY: clean


clean:
	$(RM) $(OBJ) $(EXE) $(PLTXT) $(IMG) $(TXT)
