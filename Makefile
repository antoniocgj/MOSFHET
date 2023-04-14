MOSFHET_DIR = ./
include Makefile.def

all: mosfhet

static: lib objects
	ar rcs lib/libmosfhet.a *.o && rm *.o

objects: $(SRC_MOSFHET)
	$(CC) -g -c $^ $(LIB_FLAGS) $(LIBS)

test: test/test
	./test/test

bench: test/benchmark
	./test/benchmark

mosfhet: lib lib/mosfhet

lib/mosfhet: $(SRC_MOSFHET)
	$(CC) -g -fPIC -shared -o lib/libmosfhet.so $^ $(LIB_FLAGS) $(LIBS)

lib:
	mkdir -p lib

test_debug: override TEST_FLAGS = $(DEBUG_FLAGS) -Dinline=""
test_debug: test/test
	gdb ./test/test

test/test: $(SRC_MOSFHET) test/unity_test/unity.c test/tests.c 
	$(CC) -g -o test/test $^ $(TEST_FLAGS) $(LIBS)

test/benchmark: $(SRC_MOSFHET) test/benchmark.c
	$(CC) -g -o test/benchmark $^ $(OPT_FLAGS) $(LIBS)

test_fft: test/test_fft
	./test/test_fft

test/test_fft: $(SRC_MOSFHET) test/benchmark_arith.c
	$(CC) -g -o test/test_fft $^ $(TEST_FLAGS) $(LIBS) 

clean: 
	rm --f test/test test/benchmark test/test_fft lib/libmosfhet.so lib/libmosfhet.a
