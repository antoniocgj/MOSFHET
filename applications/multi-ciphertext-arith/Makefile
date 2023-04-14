MOSFHET_DIR = ../../
include $(MOSFHET_DIR)/Makefile.def
INCLUDE_DIRS += ./include/

SRC = ufhe.c integer.c lut.c io.c ml.c
SRC_FILES = $(addprefix src/, $(SRC))

all: ufhe

test: test/test 
	$(LD_LIBS) ./test/test

bench: test/benchmark
	$(LD_LIBS) ./test/benchmark

ufhe: lib lib/ufhe

lib:
	mkdir -p lib

lib/ufhe: $(SRC_MOSFHET) $(SRC_FILES)
	$(CC) -g -fPIC -shared -o lib/libufhe.so $^ $(OPT_FLAGS) $(LIBS)

test_debug: override TEST_FLAGS = $(DEBUG_FLAGS) 
test_debug: test/test
	$(LD_LIBS) gdb ./test/test

test/test: $(SRC_MOSFHET) $(SRC_FILES) test/unity_test/unity.c test/tests.c 
	$(CC) -g -o test/test $^ $(TEST_FLAGS) $(LIBS)

test/benchmark: $(SRC_MOSFHET) $(SRC_FILES) test/benchmark.c
	$(CC) -g -o test/benchmark $^ $(OPT_FLAGS) $(LIBS)

clean: 
	rm --f test/test test/benchmark test/test_fft lib/libufhe.so
