.PHONY: default test

CFLAGS = -g -Wall -Wextra
#CFLAGS := $(CFLAGS) -DDEBUG
LFLAGS =

default: expandconfig

obj/%.o:  src/%.c
	$(CC) -c $(CFLAGS) $< -o $@

expandconfig: obj/expandconfig.o
	$(CC) $^ -o $@

clean:
	-rm -f obj/*
	-rm -f test/test.cfg.*

cleanup: clean
	-rm -f expandconfig

test: expandconfig
	./expandconfig test/test.cfg
	@cat test/test.cfg.* > test/test_result
	@if diff test/test_result test/test_sample_result > /dev/null; then echo "Test passed"; else echo "Test failed"; fi
