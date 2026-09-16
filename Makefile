
HCC=hcc

.PHONY: all clean


bin/% : %.HC
	@mkdir -p bin
	$(HCC) -o $@ $<

all: bin/ising bin/random_test

clean: 
	rm ising random_test
