
fbench: fbench.sim
	cim fbench.sim

clean:
	rm -f fbench.o fbench results.txt
	
time:	fbench
	time -p ./fbench >results.txt
	diff results.txt correct_results.txt
	
