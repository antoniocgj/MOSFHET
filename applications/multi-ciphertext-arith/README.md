# MOSFHET Multi-Ciphertext Arithmetic

This code is an example on how to use [MOSFHET](https://github.com/antoniocgj/MOSFHET) to implement the evaluation of high-level functions. 

It contains functions for arithmetic and LUT evaluation over messages decomposed and encrypted in mutiple ciphertexts. 

## Build and Running

`make` to build the shared library, `make test` to run the unit tests. You can change compiling options (e.g. change the FFT library) for MOSFHET in the [Makefile](Makefile#L28).

## Outdated techniques

Please notice that the techniques employed for arithmetic in this library are a bit outdated. Instead of using multi-ciphertext additions and multiplication, it is more efficient to:
1. Perform arithmetic over single-ciphertext messages (using native addition and tensor product with relinearization).
2. Decompose the messages using the Improved Programmable Bootstrap (PBS) to evaluate nonlinear functions with high-precision.
3. Compose the messages back with linear error growth using the multi-value extract scaling.  

MOSFHET provides optimized implementations of all these techniques (tensor product, ciphertex decomposition, improved PBS, and multi-value extract scaling). We plan to build a new library for high-level functions using them. 