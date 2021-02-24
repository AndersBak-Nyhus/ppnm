CFLAGS = -Wall -Werror -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: default
default: out.exA.txt

##out.exA.png: out.exA.txt Makefile
##	echo '\
##	set terminal png;\
##	set output "$@";\
##	set key top left;\
##	set tics out;\
##	set xlabel "x";\
##	set ylabel "y";\
##	set title "Error function";\
##	plot \
##	 "$<" using 1:2 with line title "ln(x)/sqrt(x)"\
##	' |  gnuplot

out.exA.txt: exA Makefile
	./$< > $@

.PHONEY: clean
clean:
	$(RM) out* exA
