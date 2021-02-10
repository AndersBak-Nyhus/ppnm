CC=gcc
CFLAGS=-O -std=gnu11 -march=native
LDLIBS=-lm

default: out.txt
	cat out.txt
out.txt: math
	./math  > out.txt
math: math.o
	$(CC) -o math math.o $(LDLIBS)
math.o:
	$(CC) $(CFLAGS) -c math.c
clean:
	rm -f math.o math out.txt
test:
	echo $(CC)
	echo $(CFLAGS)

