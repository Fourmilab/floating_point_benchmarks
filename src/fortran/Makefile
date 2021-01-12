
fbench:	fbench.f
	f77 -O3 fbench.f -o fbench
	
clean:
	rm -f fbench
	
time:	fbench
	echo -e "\n" | time -p ./fbench | \
		grep -v Press | grep -v Stop>results.txt
	diff results.txt correct_results.txt
