
JULIA = julia

OPTIONS = -O2 --check-bounds=no

#ITERATIONS = 110327622
ITERATIONS = 5000000

time:
	time -p $(JULIA) $(OPTIONS) fbench.jl $(ITERATIONS)

time_static:
	time -p $(JULIA) $(OPTIONS) fbench_static.jl $(ITERATIONS)
