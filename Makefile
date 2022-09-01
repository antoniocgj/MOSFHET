include Makefile.def
INCLUDE_FLAGS = $(addprefix -I, $(INCLUDE_DIRS))
FLAGS += $(INCLUDE_FLAGS)

all: mosfhet

static: lib objects
	ar rcs lib/libmosfhet.a *.o && rm *.o

objects: $(addprefix src/, $(SRC))
	$(CC) -g -c $^ $(LIB_FLAGS) $(LIBS)

test: test/test
	./test/test

bench: test/benchmark
	./test/benchmark

mosfhet: lib lib/mosfhet

lib/mosfhet: $(addprefix src/, $(SRC))
	$(CC) -g -fPIC -shared -o lib/libmosfhet.so $^ $(LIB_FLAGS) $(LIBS)

lib:
	mkdir -p lib

test_debug: override TEST_FLAGS = $(DEBUG_FLAGS) -Dinline=""
test_debug: test/test
	gdb ./test/test

test/test: $(addprefix src/, $(SRC)) test/unity_test/unity.c test/tests.c 
	$(CC) -g -o test/test $^ $(TEST_FLAGS) $(LIBS)

test/benchmark: $(addprefix src/, $(SRC)) test/benchmark.c
	$(CC) -g -o test/benchmark $^ $(OPT_FLAGS) $(LIBS)

test_fft: test/test_fft
	./test/test_fft

test/test_fft: $(addprefix src/, $(SRC)) test/benchmark_arith.c
	$(CC) -g -o test/test_fft $^ $(TEST_FLAGS) $(LIBS) 

clean: 
	rm --f test/test test/benchmark test/test_fft lib/libmosfhet.so lib/libmosfhet.a
