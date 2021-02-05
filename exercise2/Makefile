CC=gcc
CFLAGS=-O -std=gnu11 -march=native
LDLIBS=-lm

default: out.txt
	cat out.txt
out.txt: hello
	./hello > out.txt
hello: hello.o
	$(CC) -o hello hello.o $(LDLIBS)
hello.o:
	$(CC) $(CFLAGS) -c hello.c
clean:
	rm -f hello.o hello out.txt
test:
	echo $(CC)
	echo $(CFLAGS)

