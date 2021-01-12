
fbench.sh:	fbench.a68
	a68g --strict -O3 --compile fbench.a68
	
test:
	a68g fbench.a68 >results.txt
	diff results.txt ../smalltalk/results.txt
	
clean:
	rm -f fbench.sh results.txt
	
