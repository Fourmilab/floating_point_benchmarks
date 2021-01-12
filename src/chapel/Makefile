
CHPLOPTS = --fast

ITERATIONS = 10000000

TASKS = 8

serial:
	chpl $(CHPLOPTS) fbench.chpl --set partrace=false

pariter:
	chpl $(CHPLOPTS) fbench.chpl --set pariter=$(TASKS)

partrace:
	chpl $(CHPLOPTS) fbench.chpl --set partrace=true

time:
	time -p ./fbench --iterations=$(ITERATIONS)

clean:
	rm -f fbench
