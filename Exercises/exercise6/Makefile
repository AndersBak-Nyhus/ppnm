#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


erf.gnuplot.png: data.txt Makefile tabulated.erf.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated.erf.data.txt" using 1:2 with points title "tabulated data"\
	' |  gnuplot

#gamma.gnuplot.png: data.txt Makefile
#	echo '\
#	set terminal png;\
#	set output "$@";\
#	set key top left;\
#	set tics out;\
#	set xlabel "x";\
#	set ylabel "y";\
#	set title "gamma function";\
#	plot \
#	"$<" using 1:5 with line title "gamma from math.h"\
#	,"$<" using 1:6 with line title "gsl-sf-gamma"\
#	' |  gnuplot

data.txt: main
	./$< > $@

main: main.o Erf.o


.PHONEY:clean
clean:
	$(RM) *.o main *.png data.txt test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)

