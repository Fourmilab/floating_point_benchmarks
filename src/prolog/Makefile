
#ITERATIONS = 13725893

ITERATIONS = 450000

#       This Makefile is able to build three versions of the
#       benchmark:
#
#           fbench      SWI-Prolog with built-in accuracy test
#           fbench_iso  ISO Prolog with external accuracy test,
#                       built with SWI-Prolog
#           fbench_gp   ISO prolog with external accuracy test,
#                       built with GNU Prolog
#
#       All are built from the same fbench.pl source code.

all:    fbench fbench_iso fbench_gp

#       Build and time SWI-Prolog self-testing version

fbench: fbench.pl
	swipl -O -nodebug -g "run_benchmark($(ITERATIONS))." -g "halt." -o fbench -c fbench.pl

time:   FORCE
	@echo $(ITERATIONS) iterations:
	time -p ./fbench

#       Build and time ISO Prolog version, testing output with diff

fbench_iso: fbench.pl
	swipl -O -nodebug -g "run_benchmark_ISO($(ITERATIONS))." -g "halt." -o fbench_iso -c fbench.pl

time_iso:  FORCE
	@echo $(ITERATIONS) iterations:
	time -p ./fbench_iso >results.txt
	@-diff results.txt correct_results.txt && \
	    ([ $$? -eq 0 ]) || \
	    echo 'Error(s) detected in results.  This is VERY SERIOUS.'

#       Build and time GNU Prolog, testing output with diff.
#       Due to lack of tail call optimisation, the maximum
#       iteration count is limited.

fbench_gp:      fbench.pl
	gplc fbench.pl -o fbench_gp

time_gp:        FORCE
	time -p sh -c "echo 'nl, run_benchmark_ISO($(ITERATIONS)). halt.' | \
	    ./fbench_gp" | \
	    tail --lines=+6 | head -8 >results.txt
	@-diff results.txt correct_results.txt && \
	    ([ $$? -eq 0 ]) || \
	    echo 'Error(s) detected in results.  This is VERY SERIOUS.'


clean:  FORCE
	rm -f fbench fbench_iso fbench_gp results.txt

FORCE:
