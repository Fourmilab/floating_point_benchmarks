
#   Options for compiling and linking C programs

COPTS = -Wall -O3 -lm

PROGRAMS = fbench

all:	$(PROGRAMS)

fbench: fbench.c
	gcc fbench.c -o fbench $(COPTS)

clean:
	rm -f $(PROGRAMS) *.o
