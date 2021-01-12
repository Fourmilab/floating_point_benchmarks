
duh:
	@echo "What'll it be--clean, time, or timejit?"

clean:
	rm -f results.txt
	
time:
	time -p lua fbench.lua >results.txt
	diff results.txt correct_results.txt
	
timejit:
	time -p luajit fbench.lua >results.txt
	diff results.txt correct_results.txt
