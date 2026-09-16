
HCC=hcc

.PHONY: all clean

% : %.HC
	$(HCC) -o $@ $<

all: ising random_test

clean: 
	rm ising random_test
