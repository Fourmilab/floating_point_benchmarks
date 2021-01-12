
ITERATIONS = 126612394

fbench:
	scalac -opt:l:inline -opt-inline-from:** \
	       -g:none -deprecation -Xlint fbench.scala
	       
clean:
	rm -f *class

time:
	time -p scala fbench $(ITERATIONS)
	
#	Time in interpretive execution environment

timeREPL:
	time -p scala fbench.scala $(ITERATIONS)
