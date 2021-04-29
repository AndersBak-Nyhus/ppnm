CFLAGS = -Wall -Werror -std=gnu1x  -O
LDLIBS = -lm

all: out.png

out.png: out.data.txt Makefile
	echo "\
	set term png;\
	set out '$@';\
	set key top left;\
	set tics out;\
	set xlabel 'x';\
	set ylabel 'y';\
	set title 'Linear interpolation';\
	plot '$<' index 0 with points title 'data'\
	,'$<' index 1  using 1:2 with line title 'Linear interpolation'\
	" | tee log.txt |  gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt

out.data.txt: main
	./$< >out.data.txt
main: main.o linterp.o

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg test* log*
