# MOSFHET: Optimized Software for FHE over the Torus

MOSFHET is a pure-C highly-optimized implementation of [TFHE](https://github.com/tfhe/tfhe/). It includes the main techniques proposed so far for improving performance or error rate in TFHE. The library is fully portable with optional optimizations for Intel AVX2, FMA, and AVX-512. 

## Implemented Techniques

- The Functional [[1]](https://link.springer.com/chapter/10.1007/978-3-030-20951-3_20) or Programmable [[2]](https://link.springer.com/chapter/10.1007/978-3-030-78086-9_1) Bootstrap and its improved version [[3]](https://link.springer.com/chapter/10.1007/978-3-030-92078-4_23).

- The Circuit Bootstrap [[4]](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_14) and its optimizations [[3]](https://link.springer.com/chapter/10.1007/978-3-030-92078-4_23).

- The multi-value bootstrap [[3](https://link.springer.com/chapter/10.1007/978-3-030-92078-4_23),[5](https://link.springer.com/chapter/10.1007/978-3-030-12612-4_6)] and its optimizations [[6]](https://tches.iacr.org/index.php/TCHES/article/view/8793).

- The Key Switching [[7]](https://link.springer.com/chapter/10.1007/978-3-662-53887-6_1) and its optimizations [[8]](https://link.springer.com/chapter/10.1007/978-3-030-78372-3_18).

- The BlindRotate Unfolding [[9]](https://ieeexplore.ieee.org/document/8449914) and its optimizations [[10]](https://link.springer.com/chapter/10.1007/978-3-319-96878-0_17).

- The Full TRGSW bootstrap. 

- Three different approaches [[3](https://link.springer.com/chapter/10.1007/978-3-030-92078-4_23),[12](https://ia.cr/2021/1135),[13](https://ia.cr/2021/1347)] for evaluating the Full-Domain Functional Bootstrap (FDFB). 

- Public Key compression using randomness seed [[14]](10.1007/s00145-019-09319-x). 

- BFV-like multiplication [[3]](https://link.springer.com/chapter/10.1007/978-3-030-92078-4_23).

For more details, see our [paper](https://eprint.iacr.org/2022/515). 

## Build

By default, we use the AVX-512 version of SPQLIOS for fast polynomial arithmetic. It requires [AVX-512 support](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX-512). You can use the option `FFT_LIB` to specify other libraries. 

Default compilation (using AVX-512 SPQLIOS):

> `make`

SPQLIOS (FMA):

> ```make FFT_LIB=spqlios```

FFNT library (pure-C, fully portable):

> ```make FFT_LIB=ffnt```

For other compiling options, see the [Makefile](Makefile). 

## Running

There are two main ways of using MOSFHET:
1. The most efficient is to compile your code and MOSFHET together. We do that for our [benchmark.c](test/benchmark.c) and [tests.c](test/tests.c) files. See the `test/benchmark` rule in the [Makefile](Makefile#L48). 
2. Dynamic Link. After compiling MOSFHET as a shared library, you can dynamically link it with your code. See [MOSFHET_MCA](https://github.com/antoniocgj/MOSFHET_MCA) for an example. 

## Examples

For examples on how to use MOSFHET, see our [unit tests file](test/tests.c) and [MOSFHET_MCA](https://github.com/antoniocgj/MOSFHET_MCA). 

## Unit tests and Benchmark

We provide a set of [unit tests](test/tests.c) and a simple [benchmark](test/benchmark.c) file for the library. They use parameters hard-coded at the beginning of each file. The default parameters are high memory consuming (they are the same as TFHEpp Level 2). The parameters can be reduced in exchange for performance (especially for the Key Switching) or error rate.

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

```
@misc{cryptoeprint:2022:515,
    author       = {Antonio Guimarães and
		    Edson Borin and
		    Diego F. Aranha},
    title        = {MOSFHET: Optimized Software for FHE over the Torus},
    howpublished = {Cryptology ePrint Archive, Report 2022/515},
    year         = {2022},
    note         = {\url{https://ia.cr/2022/515}},
}
```

[The paper](https://eprint.iacr.org/2022/515) consider the [initial commit (0d58320559)]( https://github.com/antoniocgj/MOSFHET/tree/0d5832055900d1376f7dadbbf5093b911e96a7fc) of the library in this repository. 

## License

[Apache License Version 2.0](LICENSE)

This repository includes code from the following third party libraries:
- [FFNT](https://gitlab.fit.cvut.cz/klemsjak/ffnt-benchmark): [MIT License](https://gitlab.fit.cvut.cz/klemsjak/ffnt-benchmark/blob/master/LICENSE), Copyright (c) 2021 Jakub Klemsa
- [SPQLIOS](https://github.com/tfhe/tfhe/tree/master/src/libtfhe/fft_processors/spqlios): [Apache License Version 2.0](https://github.com/tfhe/tfhe/blob/master/LICENSE), Copyright 2016 - Nicolas Gama <nicolas.gama@gmail.com> et al.
- [Unity Test](https://github.com/ThrowTheSwitch/Unity): [MIT License](https://github.com/ThrowTheSwitch/Unity/blob/master/LICENSE.txt), Copyright (c) <year> 2007-21 Mike Karlesky, Mark VanderVoord, Greg Williams
- [FIPS202 from Kyber](https://github.com/pq-crystals/kyber/blob/master/ref/fips202.c): [Public Domain](https://creativecommons.org/share-your-work/public-domain/cc0/)
- [xoshiro / xoroshiro](https://prng.di.unimi.it/): [Public Domain](https://creativecommons.org/share-your-work/public-domain/cc0/), David Blackman and Sebastiano Vigna (vigna@acm.org)


Additionally, our library may contain small code snippets, variable names, or implementation logic based on or adapted from:
- [TFHE](https://github.com/tfhe/tfhe/): [Apache License Version 2.0](https://github.com/tfhe/tfhe/blob/master/LICENSE), Copyright 2016 - Nicolas Gama <nicolas.gama@gmail.com> et al.
- [TFHEpp](https://github.com/virtualsecureplatform/TFHEpp): [Apache License Version 2.0](https://github.com/virtualsecureplatform/TFHEpp/blob/master/LICENSE), Copyright 2019 Kotaro MATSUOKA
- [dbush](https://stackoverflow.com/questions/48043811/creating-a-function-to-check-if-malloc-succeeded): [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0)
- [Lattigo](https://github.com/tuneinsight/lattigo): [Apache License Version 2.0](https://github.com/tuneinsight/lattigo/blob/master/LICENSE)
