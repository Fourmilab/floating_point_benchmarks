
FBC = fbc
#FBC = $(HOME)/linuxtools/FreeBASIC-1.05.0-linux-x86_64/bin/fbc

FBOPTS =

ITERATIONS = 1000000
#ITERATIONS = 127172531

fbench: fbench.bas
	$(FBC) $(FBOPTS) fbench.bas

time:   fbench
	@echo $(ITERATIONS) iterations:
	time -p ./fbench $(ITERATIONS) >results.txt
	@-diff results.txt correct_results.txt && \
	    ([ $$? -eq 0 ]) || \
	    echo 'Error(s) detected in results.  This is VERY SERIOUS.'

clean:
	rm -f fbench results.txt

FORCE:
