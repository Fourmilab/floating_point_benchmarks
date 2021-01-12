
GM2 = gm2 -O3

LIBS = -fiso -flibs=iso,pim -Wpedantic

ITERATIONS = 10000000
#ITERATIONS = 175495723

fbench: fbench.mod
	$(GM2) $(LIBS) fbench.mod -o fbench

time:   fbench
	@echo $(ITERATIONS) iterations:
	time -p ./fbench $(ITERATIONS) >results.txt
	@-diff results.txt correct_results.txt && \
	    ([ $$? -eq 0 ]) || \
	    echo 'Error(s) detected in results.  This is VERY SERIOUS.'

clean:
	rm -f *.o *.s *.cpp *.lst fbench results.txt
