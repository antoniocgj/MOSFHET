# MOSFHET: Optimized Software for FHE over the Torus

MOSFHET is a pure-C highly-optimized implementation of [TFHE](https://github.com/tfhe/tfhe/). It includes the main techniques proposed so far for improving performance or error rate in TFHE. The library is fully portable with optional optimizations for Intel AVX2, FMA, and AVX-512. 

<!-- ## Implemented Techniques

- The Functional [[1]](#1) or Programmable [[2]](#1) Bootstrap and its improved version [[3]](#1).

- The Circuit Bootstrap [[4]](#1) and its optimizations [[3]](#1).

- The multi-value bootstrap [[3,5]](#5) and its optimizations [[6]](#1).

- The Key Switching [[7]](#1) and its optimizations [[8]](#1).

- The BlindRotate Unfolding [[9]](#9) and its optimizations [[10]](#1).

- The Full TRGSW bootstrap [[11]](#1). 

- Three different approaches [[3,12,13]](#1) for evaluating the Full-Domain Functional Bootstrap (FDFB). 

- Public Key compression using randomness seed [[14]](#1). 

- BFV-like multiplication [[3]](#1). -->

## Build and Running

By default, we use the AVX-512 version of SPQLIOS for fast polynomial arithmetic. It requires [AVX-512 support](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX-512). You can use the option `FFT_LIB` to specify other libraries. 

Default compilation (using AVX-512 SPQLIOS):

> `make`

SPQLIOS (FMA):

> ```make FFT_LIB=spqlios```

FFNT library (pure-C, fully portable):

> ```make FFT_LIB=ffnt```

For other compiling options, see the [Makefile](Makefile). 

### Unit tests and Benchmark

We provide a set of [unit tests](test/tests.c) and a simple [benchmark](test/benchmark.c) file for the library. They use parameters hard-coded at the begining of each file. The default parameters are high memory consuming (they are the same as TFHEpp Level 2). The parameters can be reduced in exchange for performance (especially for the Key Switching) or error rate.

To run the unit tests:

> `make test -B`
> 
> `make test FFT_LIB=spqlios -B`
> 
> `make test FFT_LIB=ffnt -B` 

To run the benchmark:

> `make bench -B`
> 
> `make bench FFT_LIB=spqlios -B`
> 
> `make bench FFT_LIB=ffnt -B` 

## Citation

tbd

## License

[Apache License Version 2.0](LICENSE)

This repository includes code from the following third party libraries:
- [FFNT](https://gitlab.fit.cvut.cz/klemsjak/ffnt-benchmark): [MIT License](https://gitlab.fit.cvut.cz/klemsjak/ffnt-benchmark/blob/master/LICENSE), Copyright (c) 2021 Jakub Klemsa
- [SPQLIOS](https://github.com/tfhe/tfhe/tree/master/src/libtfhe/fft_processors/spqlios): [Apache License Version 2.0](https://github.com/tfhe/tfhe/blob/master/LICENSE), Copyright 2016 - Nicolas Gama <nicolas.gama@gmail.com> et al.
- [Unity Test](https://github.com/ThrowTheSwitch/Unity): [MIT License](https://github.com/ThrowTheSwitch/Unity/blob/master/LICENSE.txt), Copyright (c) <year> 2007-21 Mike Karlesky, Mark VanderVoord, Greg Williams
- [FIPS202 from Kyber](https://github.com/pq-crystals/kyber/blob/master/ref/fips202.c): [Public Domain](https://creativecommons.org/share-your-work/public-domain/cc0/)
- [xoshiro / xoroshiro](https://prng.di.unimi.it/): [Public Domain](https://creativecommons.org/share-your-work/public-domain/cc0/), David Blackman and Sebastiano Vigna (vigna@acm.org)


Additionally, our library may contain small snippets of code, variable names, or implementation logic based on or adapted from:
- [TFHE](https://github.com/tfhe/tfhe/): [Apache License Version 2.0](https://github.com/tfhe/tfhe/blob/master/LICENSE), Copyright 2016 - Nicolas Gama <nicolas.gama@gmail.com> et al.
- [TFHEpp](https://github.com/virtualsecureplatform/TFHEpp): [Apache License Version 2.0](https://github.com/virtualsecureplatform/TFHEpp/blob/master/LICENSE), Copyright 2019 Kotaro MATSUOKA
- [dbush](https://stackoverflow.com/questions/48043811/creating-a-function-to-check-if-malloc-succeeded): [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0)
- [Lattigo](https://github.com/tuneinsight/lattigo): [Apache License Version 2.0](https://github.com/tuneinsight/lattigo/blob/master/LICENSE)