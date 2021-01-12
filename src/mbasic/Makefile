
#       Makefile to run with Bas:

#       http://www.moria.de/~michael/bas/

BAS = bas
#BAS = $(HOME)/linuxtools/bas-2.4/bas

#ITERATIONS = 3056858
ITERATIONS = 100000

time:
	echo "$(ITERATIONS)\\n\\n" | \
		sh -c "time -p $(BAS) fbench.bas" >results.txt
	@-diff results.txt correct_results.txt && \
	    ([ $$? -eq 0 ]) || \
	    echo 'Error(s) detected in results.  This is VERY SERIOUS.'
	rm results.txt

