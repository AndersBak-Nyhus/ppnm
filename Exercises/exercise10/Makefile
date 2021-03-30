CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex
	pdflatex $<
	pdflatex $<

out.data.txt: main
	./$< > $@

main: main.o ex.o

fig-gpl.tex: fig-gpl.gpi out.data.txt
	gnuplot $<


clean:
	$(RM) ms.pdf *.log *.aux out* main *.o 
