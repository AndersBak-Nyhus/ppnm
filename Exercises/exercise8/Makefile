CFLAGS = -Wall -pthread -Ofast -march=native


Out.error.png: Out.data.txt Makefile
	echo '\
	set terminal png ;\
	set output "$@" ;\
	set title "Monte Carlo calculation of $$\pi$$" ;\
	set logscale x;\
	set logscale y;\
	set xlabel "number of points $$N$$";\
	set ylabel "absolute error";\
	a=1;\
	f(x)=a/sqrt(x);\
	fit f(x) "$<" using 1:3 via a;\
	plot [][1e-5:8e-3]\
	"$<" using 1:3 with linespoints title "actual error"\
	,f(x) with lines title "$$1/sqrt{N}$$"\
	' | gnuplot

ns = 1 1.5 2 4 8

Out.data.txt: main Makefile
	>$@
	for n in $(ns); do ./$< "$$n"e5 >> $@; done
	for n in $(ns); do ./$< "$$n"e6 >> $@; done
	for n in $(ns); do ./$< "$$n"e7 >> $@; done


.PHONEY: clean

clean:
	$(RM) main Out* *.log
