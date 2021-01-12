
CPP = g++

#COPTS = -g -Wall
COPTS = -O3 -Wall

#   Iterations to run with time target
ITERATIONS = 1000000

#   Precision in bits to use in fbench_mpfr builds
MPFR_PRECISION = 128

PROGRAMS = fbench fbench_ld fbench_128 fbench_mpfr

#   Standard version, using "double"

fbench: fbench.cpp
	$(CPP) $(COPTS) fbench.cpp -o fbench -lm

#   Version using "long double"

fbench_ld: fbench.cpp
	$(CPP) $(COPTS) -DLONG_DOUBLE=1 fbench.cpp -o fbench_ld -lm

#   Version using GCC's libquadmath 128-bit floating point

fbench_128: fbench.cpp
	$(CPP) $(COPTS) -DFLOAT128=1 fbench.cpp -o fbench_128 -lquadmath

#   Version using the MPFR multiple precision package

fbench_mpfr: fbench.cpp
	$(CPP) $(COPTS) -DFLOAT_MPFR=1 -DMPFR_PRECISION=$(MPFR_PRECISION) \
                fbench.cpp -o fbench_mpfr -lgmp -lmpfr

all:    $(PROGRAMS)

clean:
	rm -f $(PROGRAMS) core*

time:   fbench
	time -p ./fbench $(ITERATIONS)
