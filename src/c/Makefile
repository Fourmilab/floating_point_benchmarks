
#   Options for compiling and linking C fbench programs

#       For GCC
CC = gcc -Wall
#       For CLANG
#CC = $(HOME)/linuxtools/clang+llvm-5.0.0-linux-x86_64-ubuntu16.04/bin/clang

#       C compiler common options
COPTS = -O3 -lm

#       Enable to use internal trigonometric functions in fbench
#INTRIG = -DINTRIG

#       Iterations to run in timing tests
ITERATIONS_FBENCH = 1000000
ITERATIONS_FFBENCH = 1000

PROGRAMS = fbench ffbench fbench_ansi ffbench_ansi

all:	$(PROGRAMS)

fbench: fbench.c
	$(CC) $(INTRIG) fbench.c -o fbench $(COPTS)

fbench_ansi: fbench_ansi.c
	$(CC) $(INTRIG) fbench_ansi.c -o fbench_ansi $(COPTS)

ffbench: ffbench.c
	$(CC) ffbench.c -o ffbench $(COPTS)

ffbench_ansi: ffbench_ansi.c
	$(CC) ffbench_ansi.c -o ffbench_ansi $(COPTS)

clean:
	rm -f $(PROGRAMS) *.o *.out *.jcp core*

time_fbench: fbench
	echo -e "" | sh -c "time -p ./fbench $(ITERATIONS_FBENCH)"

time_fbench_ansi: fbench_ansi
	echo -e "" | sh -c "time -p ./fbench_ansi $(ITERATIONS_FBENCH)"

time_ffbench: ffbench
	time -p ./ffbench $(ITERATIONS_FFBENCH)

time_ffbench_ansi: ffbench_ansi
	time -p ./ffbench_ansi $(ITERATIONS_FFBENCH)
