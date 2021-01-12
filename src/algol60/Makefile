
time:
	marst fbench.alg >fbench.c
	gcc fbench.c -O3 --static -o fbench -lalgol -lm
	time -p ./fbench >output.txt
	perl reformat.pl output.txt
	
clean:
	rm -f fbench.c fbench output.txt
