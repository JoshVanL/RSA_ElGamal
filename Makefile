help:
	# all:         run all
	# build:       build file
	# build_debug: build file in debug mode (ElGamal Encrypt uses k=1)
	# test:        run tests
	# check:       build in debug mode then test
	# clean:       clean build

all: build check clean

build:
	gcc -Wall -std=gnu99 -O3 -o modmul -lgmp modmul.c

build_debug:
	gcc -Wall -std=gnu99 -O3 -o modmul -lgmp -DDEBUG modmul.c

check: build_debug test

test:
	for number in 1 2 3 4 ; do \
		./modmul stage$$number < stage$$number.input > foo; \
		cmp stage$$number.output foo || echo "stage $$number failed"; \
	done
	rm -f foo


clean:
	rm -f core modmul foo
