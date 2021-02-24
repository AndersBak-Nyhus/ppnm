CC = clang
CFLAGS = -Wall

.PHONEY: default
default: out.txt
	cat $<

out.txt: main
	./$< > $@

.PHONEY: clean
clean:
	$(RM) main out.txt
