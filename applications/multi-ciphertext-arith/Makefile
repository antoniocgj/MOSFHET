MAKEFILE := $(shell git rev-parse --show-toplevel)/Makefile.def
include $(MAKEFILE)
INCLUDE_FLAGS = $(addprefix -I../../, $(INCLUDE_DIRS)) -I./include
FLAGS += $(INCLUDE_FLAGS)
MOSFHET_SRC = $(addprefix ../../src/, $(SRC))

MY_SRC = ufhe.c integer.c lut.c io.c ml.c
SRC_FILES = $(addprefix src/, $(MY_SRC))
SRC_FILES += $(MOSFHET_SRC)

all: main

main:  main.c
	$(CC) -g -o main $^ $(OPT_FLAGS) $(LIBS)

all: ufhe

test: test/test 
	$(LD_LIBS) ./test/test

bench: test/benchmark
	$(LD_LIBS) ./test/benchmark

ufhe: lib lib/ufhe

lib:
	mkdir -p lib

lib/ufhe: $(SRC_FILES)
	$(CC) -g -fPIC -shared -o lib/libufhe.so $^ $(OPT_FLAGS) $(LIBS)

test_debug: override TEST_FLAGS = $(DEBUG_FLAGS) 
test_debug: test/test
	$(LD_LIBS) gdb ./test/test

test/test: $(SRC_FILES) test/unity_test/unity.c test/tests.c 
	$(CC) -g -o test/test $^ $(TEST_FLAGS) $(LIBS)

test/benchmark: $(SRC_FILES) test/benchmark.c
	$(CC) -g -o test/benchmark $^ $(OPT_FLAGS) $(LIBS)

clean: 
	rm --f test/test test/benchmark test/test_fft lib/libufhe.so
