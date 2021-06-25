CFLAGS = -Wall
LDLIBS = -lm

.PHONEY: all clean
all              : out.txt           ; cat $<
out.txt          : main              ; ./$< > $@
main             : main.o komplex.o equal.o 
main.o komplex.o : komplex.h        
equal.o					 : equal.c funcs.h
clean            :                   ; $(RM) main main.o komplex.o out.txt
