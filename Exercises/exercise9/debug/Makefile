CFLAGS = -O -w -std=gnu99
CFLAGS += $(gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
LDLIBS += -lm -lgsl

all: output.txt ; cat $<
output.txt: main ; ./$< > $@
clean: ; $(RM) main output*
