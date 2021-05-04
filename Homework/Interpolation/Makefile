CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

TXT = output.txt
LIN = Lin_output.txt

all: out.png

out.png: $(TXT) $(LIN) Makefile
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
	,'$<' using 1:3 with line title 'integral'\
	,'$<' using 1:4 with line title 'gsl Linear interpolation'\
	,'$<' using 1:5 with line title 'gsl integration'\
	" | tee log.txt |  gnuplot
	sed -i 's/;/;\n/g;s/,/,\n/g' log.txt


$(LIN)  $(TXT): main
	./$< $(LIN) 1>$(TXT)

main: main.o linterp.o

.PHONEY: clean
clean:
	$(RM) *.o out* *.svg test* log*
