help:
	# all:   build and run tests
	# build: build file
	# check: run tests
	# clean: clean build

all: build check

build:
	gcc -Wall -std=gnu99 -O3 -o modmul -lgmp -lm modmul.c

check:
	for number in 1 2 3 4 ; do \
		./modmul stage$$number < stage$$number.input > foo; \
		cmp stage$$number.output foo || echo "stage $$number failed"; \
	done
	rm -f foo


clean:
	rm -f core modmul foo
