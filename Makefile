CC = gcc-10
FLAGS = -Wall -Wno-unused-function -Wno-unused-result -funroll-all-loops -march=native -lm -I./include/
DEBUG_FLAGS = -O0 -g $(FLAGS)
OPT_FLAGS = -O3 -fwhole-program -flto -DNDEBUG $(FLAGS)
LIB_FLAGS = -O3 -DNDEBUG $(FLAGS) -flto -fuse-ld=gold -fuse-linker-plugin
LIBS = 
LD_LIBS =
BUILD_LIBS = 
TEST_FLAGS = $(OPT_FLAGS)
FFT_LIB = spqlios_avx512

SRC = keyswitch.c bootstrap.c tlwe.c trlwe.c trlwe_compressed.c trgsw.c misc.c	polynomial.c register.c sha3/fips202.c fft/karatsuba.c

ifeq ($(FFT_LIB),spqlios)
	FLAGS += -DUSE_SPQLIOS
	SRC += ./fft/spqlios/spqlios-fft-fma.s ./fft/spqlios/spqlios-ifft-fma.s ./fft/spqlios/spqlios-fft-impl.c ./fft/spqlios/fft_processor_spqlios.c
else ifeq ($(FFT_LIB),spqlios_avx512)
	FLAGS += -DUSE_SPQLIOS -DAVX512_OPT
	SRC += ./fft/spqlios/spqlios-fft-avx512.s ./fft/spqlios/spqlios-ifft-avx512.s ./fft/spqlios/spqlios-fft-impl-avx512.c ./fft/spqlios/fft_processor_spqlios.c
else
	FLAGS += -DPORTABLE_BUILD
 	SRC += ./fft/ffnt/ffnt.c
endif

all: mosfhet

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
	rm --f test/test test/benchmark test/test_fft lib/libmosfhet.so
