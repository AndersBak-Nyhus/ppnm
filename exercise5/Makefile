CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY:default
default:out.cmdline.txt

out.cmdline.txt:cmdline input.txt
	./cmdline $$(cat input.txt) > $@

input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY:clean
clean:
	$(RM) input.txt out* cmdline
