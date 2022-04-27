#ifndef FFNT_H
#define FFNT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

extern void * safe_malloc(size_t size);
void * safe_aligned_malloc(size_t size);
// #define FMA_OPT
#ifdef FMA_OPT
#include <immintrin.h>
#endif

// #define USE_LONG_DOUBLE

/*******************************************************************************
 * Used floating point type.
 *
 * double: IEEE 754 type, 52+1-bit precision, 64-bit total
 * long double: IEEE 754 type, 64-bit precision, 80-bit total
 * */
typedef
#ifdef USE_LONG_DOUBLE
long double
#else
double
#endif
fp_t;

/*******************************************************************************
 * Type of tables.
 * */
typedef enum {
    FFT_e   = 0,
    IFFT_e  = 1,
} transform_t;

/*******************************************************************************
 * Structure that holds FFT tables.
 * */
typedef struct FftTables
{
    transform_t type;
    size_t n;
    size_t cos_denom;
    size_t tgt_size;
    size_t * bit_reversed;
    fp_t * cos_table;
    fp_t * sin_table;
    fp_t * ct_path_tables;
    fp_t * gs_path_tables;
} FftTables;

/*******************************************************************************
 * Initialize FFT tables.
 * */
FftTables *  fft_init(const size_t n);
/*******************************************************************************
 * Initialize inverse FFT tables.
 * */
FftTables * ifft_init(const size_t n);

/*******************************************************************************
 * Run FFT (via Cooley-Tukey data path).
 * */
void  fft_transform(const FftTables *const fft_tables,
                    fp_t *const real,
                    fp_t *const imag);
/*******************************************************************************
 * Run inverse FFT (via Cooley-Tukey data path).
 * */
void ifft_transform(const FftTables *const ifft_tables,
                    fp_t *const real,
                    fp_t *const imag);

/*******************************************************************************
 * Run FFNT (does not apply bit-reverse).
 * */
void ffnt_transform(const FftTables *const ffnt_tables_2N,
                    const FftTables *const ffnt_tables_N_2,
                    fp_t *const real,
                    fp_t *const imag);
/*******************************************************************************
 * Run inverse FFNT (does not apply bit-reverse).
 * */
void iffnt_transform(const FftTables *const ffnt_tables_2N,
                     const FftTables *const iffnt_tables_N_2,
                     fp_t *const real,
                     fp_t *const imag);

/*******************************************************************************
 * Run FFT via Cooley-Tukey data path.
 *
 * Choose whether to use bit-reverse of the input using the +bitrev+ switch.
 *
 * */
void  fft_transform_CT(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run inverse FFT via Cooley-Tukey data path.
 *
 * Choose whether to use bit-reverse of the input using the +bitrev+ switch.
 *
 * */
void ifft_transform_CT(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run FFT via Gentleman-Sande data path.
 *
 * Choose whether to use bit-reverse of the output using the +bitrev+ switch.
 *
 * */
void  fft_transform_GS(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev);
/*******************************************************************************
 * Run inverse FFT via Gentleman-Sande data path.
 *
 * Choose whether to use bit-reverse of the output using the +bitrev+ switch.
 *
 * */
void ifft_transform_GS(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev);

/*******************************************************************************
 * Final cleanup.
 * */
void tables_destroy(FftTables *const tables);


/* Functions for TFHE */

typedef struct{
  double *fpr, *fpi;
  const FftTables * ffnt_2n_tables, * fft_n_2_tables, * ifft_n_2_tables;
  int N;
} * FFT_Processor_FFNT;

FFT_Processor_FFNT new_FFT_Processor_FFNT(int N);
void execute_reverse_torus64(double * res, const uint64_t * a, FFT_Processor_FFNT proc);
void execute_direct_torus64(uint64_t * res, const double * a, FFT_Processor_FFNT proc);

#ifdef __cplusplus
}
#endif

#endif   // FFNT_H
