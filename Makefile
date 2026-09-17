
HCC=hcc

.PHONY: all clean test


bin/% : %.HC
	@mkdir -p bin
	$(HCC) -o $@ $<


bin/test/% : test/%.HC
	@mkdir -p bin/test
	$(HCC) -o $@ $<

all: bin/ising

test: bin/test/test_rng bin/test/test_rot

clean: 
	rm -rf bin
