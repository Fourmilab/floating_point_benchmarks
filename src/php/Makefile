
ITERATIONS = 100000
#ITERATIONS = 32952075

PHP = php

duh:
	@echo "What'll it be?  time clean php_version"

time:
	$(PHP) fbench.php $(ITERATIONS) >fbench.html

clean:
	rm -f fbench.html

php_version:
	$(PHP) -v
