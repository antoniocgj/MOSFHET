#include <math.h>
#include <tgmath.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "ffnt.h"

//~ #define VERBOSE
#ifdef VERBOSE
    #include <stdio.h>
    #include <io.h>
#endif   // #ifdef VERBOSE


#ifdef __cplusplus
extern "C" {
#endif


// ===   Macros   ==============================================================

// #ifndef M_PIl
//     #define M_PIl 3.141592653589793238462643383279502884L
// #endif


// ===   Private function prototypes   =========================================

/*******************************************************************************
 *  Implements FFT tables initialization.
 */
static FftTables * omegas_init(const size_t n,
                               const bool inverse);

/*******************************************************************************
 *  Implements FFT via Cooley-Tukey data path.
 */
static void dft_CT(const FftTables *const tables,
                   fp_t *const real,
                   fp_t *const imag,
                   const bool inverse,
                   const bool bitrev);

/*******************************************************************************
 *  Implements FFT via Gentleman-Sande data path.
 */
static void dft_GS(const FftTables *const tables,
                   fp_t *const real,
                   fp_t *const imag,
                   const bool inverse,
                   const bool bitrev);

/*******************************************************************************
 *  Inline function to permute +re+ and +im+ array, given a permutation.
 */
static inline void permute_ary(const size_t *const permutation,
                               fp_t *const re,
                               fp_t *const im,
                               const size_t n,
                               const bool inverse);
static fp_t accurate_sine(uint64_t i, uint64_t n);
static int32_t floor_log2(size_t n);
static size_t reverse_bits(const size_t x, const uint32_t n);


// ===   Function implementations   ============================================

// ---   Init tables   ---------------------------------------------------------

FftTables * fft_init(const size_t n)
{
    return omegas_init(n, false);
}
FftTables * ifft_init(const size_t n)
{
    return omegas_init(n, true);
}

static FftTables * omegas_init(const size_t n,
                               const bool inverse)
{
    //TODO check these checks
    // check size
    if (n < 4 || (n & (n - 1)) != 0)
        return NULL;  // too small or not a power of 2
    if ((n / 2 + 2 * n - 8) > SIZE_MAX / sizeof(fp_t) || n > SIZE_MAX / sizeof(size_t))
        return NULL;  // too large

    // alloc table structure
    FftTables * tables = (FftTables *)safe_aligned_malloc(sizeof(FftTables));
    if (tables == NULL)
        return NULL;

    tables->n = n;
    tables->type = inverse ? IFFT_e : FFT_e;
    tables->cos_denom = n;
    tables->tgt_size = 2 * (tables->n - 4);

    // alloc arrays
    tables->bit_reversed = (size_t *)safe_aligned_malloc(tables->n * sizeof(size_t));

    // cos/sin tables (needed for twisting in FFNT)
    tables->cos_table   = (fp_t *)safe_aligned_malloc(tables->cos_denom / 2 * sizeof(fp_t));
    tables->sin_table   = (fp_t *)safe_aligned_malloc(tables->cos_denom / 2 * sizeof(fp_t));
    tables->ct_path_tables = (fp_t *)safe_aligned_malloc(tables->tgt_size      * sizeof(fp_t));
    tables->gs_path_tables = (fp_t *)safe_aligned_malloc(tables->tgt_size      * sizeof(fp_t));

    // check allocation
    if (tables->bit_reversed    == NULL ||
        tables->cos_table       == NULL ||
        tables->sin_table       == NULL ||
        tables->ct_path_tables  == NULL ||
        tables->gs_path_tables  == NULL)
    {
        free(tables->bit_reversed);
        free(tables->cos_table);
        free(tables->sin_table);
        free(tables->ct_path_tables);
        free(tables->gs_path_tables);
        free(tables);
        return NULL;
    }

    // bit-reverse permutation
    size_t i;
    int32_t levels = floor_log2(tables->n);
    for (i = 0; i < tables->n; i++)
        tables->bit_reversed[i] = reverse_bits(i, levels);

    // cos/sin tables
    for (i = 0; i < tables->cos_denom / 2; i++)
    {
        tables->cos_table[i] = accurate_sine(i + (tables->cos_denom) / 4, (tables->cos_denom));
        tables->sin_table[i] = (inverse ? -1 : 1) * accurate_sine(i, (tables->cos_denom));
    }

    size_t size;
    size_t j, k;

    // trigonometric tables for each FFT internal level of the CT data path
    k = 0;
    for (size = 8; size <= tables->n; size <<= 1)
    {
        for (i = 0; i < size / 2; i += 4)
        {
            for (j = 0; j < 4; j++, k++)
            {
                tables->ct_path_tables[k] = accurate_sine(i + j + size / 4, size);  // Cosine
            }

            for (j = 0; j < 4; j++, k++)
            {
                tables->ct_path_tables[k] = (inverse ? -1 : 1) * accurate_sine(i + j, size);  // Sine
            }
        }
        // prevent overflow
        if (size == tables->n)
            break;
    }

    // trigonometric tables for each FFT internal level of the GS data path
    k = 0;
    for (size = tables->n; size >= 8; size >>= 1)
    {
        for (i = 0; i < size / 2; i += 4)
        {
            for (j = 0; j < 4; j++, k++)
            {
                tables->gs_path_tables[k] = accurate_sine(i + j + size / 4, size);  // Cosine, TODO
            }

            for (j = 0; j < 4; j++, k++)
            {
                tables->gs_path_tables[k] = (inverse ? -1 : 1) * accurate_sine(i + j, size);  // Sine, TODO
            }
        }
    }

    return tables;
}


// ---   FFT   -----------------------------------------------------------------

// in general, use the CT data path
void  fft_transform(const FftTables *const fft_tables,
                    fp_t *const real,
                    fp_t *const imag)
{
    fft_transform_CT(fft_tables, real, imag, true);
}
void ifft_transform(const FftTables *const ifft_tables,
                    fp_t *const real,
                    fp_t *const imag)
{
    ifft_transform_CT(ifft_tables, real, imag, true);
}

// call specific data path
void  fft_transform_CT(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev)
{
    dft_CT(fft_tables, real, imag, false, bitrev);
}
void ifft_transform_CT(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev)
{
    dft_CT(fft_tables, real, imag, true, bitrev);
}
void  fft_transform_GS(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev)
{
    dft_GS(fft_tables, real, imag, false, bitrev);
}
void ifft_transform_GS(const FftTables *const fft_tables,
                       fp_t *const real,
                       fp_t *const imag,
                       const bool bitrev)
{
    dft_GS(fft_tables, real, imag, true, bitrev);
}


// ---   FFNT   ----------------------------------------------------------------

void ffnt_transform(const FftTables *const ffnt_tables_2N,
                    const FftTables *const ffnt_tables_N_2,
                    fp_t *const real,
                    fp_t *const imag)
{
    const fp_t *const fan_cos = ffnt_tables_2N->cos_table;
    const fp_t *const fan_sin = ffnt_tables_2N->sin_table;
    size_t n_2 = ffnt_tables_N_2->n;

    // fan-style multiplication
    for (size_t i = 0; i < n_2; i++)
    {
        imag[i] = real[i] * fan_sin[i] + real[n_2+i] * fan_cos[i];
        real[i] = real[i] * fan_cos[i] - real[n_2+i] * fan_sin[i];
    }

    fft_transform_GS(ffnt_tables_N_2, real, imag, false);
}

void iffnt_transform(const FftTables *const ffnt_tables_2N,
                     const FftTables *const iffnt_tables_N_2,
                     fp_t *const real,
                     fp_t *const imag)
{
    const fp_t *const fan_cos = ffnt_tables_2N->cos_table;
    const fp_t *const fan_sin = ffnt_tables_2N->sin_table;
    size_t n_2 = iffnt_tables_N_2->n;

    ifft_transform_CT(iffnt_tables_N_2, real, imag, false);

    // fan-style multiplication (negative index)
    for (size_t i = 0; i < n_2; i++)
    {
        real[n_2+i] = -real[i] * fan_sin[i] + imag[i] * fan_cos[i];
        real[i]     =  real[i] * fan_cos[i] + imag[i] * fan_sin[i];
    }
}


// ---   Implementation of FFT via the Cooley-Tukey data path   ----------------

static void dft_CT(const FftTables *const tables,
                   fp_t *const real,
                   fp_t *const imag,
                   const bool inverse,
                   const bool bitrev)
{
    // check table type
    if (!((inverse && tables->type == IFFT_e) || (!inverse && tables->type == FFT_e))) return;

    const size_t n = tables->n;

    size_t i, j;
    fp_t tpre, tpim;

    size_t halfsize, tablestep, size = 2;

    // Cooley-Tukey decimation-in-time
    const fp_t * trigtables = tables->ct_path_tables;   // no *const since the pointer is iterated

    if (bitrev)
    {
#ifdef VERBOSE
    printf("\nBefore bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        const size_t *const bitreversed = tables->bit_reversed;
        permute_ary(bitreversed, real,
                    imag,
                    n, inverse);

#ifdef VERBOSE
    printf("\nAfter bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE
    }
    else if (inverse)   // n.b., division is a part of bitreverse
    {
        for (i = 0; i < n; i++)
        {
            real[i] /= n;
            imag[i] /= n;
        }
    }

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(2));
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 2 merge (special)
    if (n >= 2)
    {
        for (i = 0; i < n; i += 2)
        {
            // addition and subtraction
            tpre = real[i];
            tpim = imag[i];

            real[i] += real[i + 1];
            imag[i] += imag[i + 1];
            real[i + 1] = tpre - real[i + 1];
            imag[i + 1] = tpim - imag[i + 1];
        }
    }

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(4));
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 4 merge (special)
    if (n >= 4)
    {
        if (!inverse)
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices
                tpre = real[i + 2];   // y
                tpim = imag[i + 2];
                real[i + 2] = real[i] - tpre;   // y = x - y
                imag[i + 2] = imag[i] - tpim;
                real[i] += tpre;   // x = x + y
                imag[i] += tpim;

                // odd indices (times -i)
                tpre =  imag[i + 3];   // -iy
                tpim = -real[i + 3];
                real[i + 3] = real[i + 1] - tpre;   // y = x - iy
                imag[i + 3] = imag[i + 1] - tpim;
                real[i + 1] += tpre;   // x = x + iy
                imag[i + 1] += tpim;
            }
        }
        else
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices
                tpre = real[i + 2];
                tpim = imag[i + 2];

                real[i + 2] = real[i] - tpre;
                imag[i + 2] = imag[i] - tpim;
                real[i] += tpre;
                imag[i] += tpim;

                // odd indices (times i)
                tpre = -imag[i + 3];
                tpim =  real[i + 3];

                real[i + 3] = real[i + 1] - tpre;
                imag[i + 3] = imag[i + 1] - tpim;
                real[i + 1] += tpre;
                imag[i + 1] += tpim;
            }
        }
    }
#ifdef FMA_OPT
    __m256d * trigtablesv = (__m256d *) trigtables;
    __m256d * realv = (__m256d *) real;
    __m256d * imagv = (__m256d *) imag;
    for (size = 2; size <= n/4; size <<= 1)
    {
        halfsize = size >> 1;

        for (i = 0; i < n/4; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 1, tablestep += 2){
                uint64_t vi = i + j;        // Vector index
                uint64_t ti = tablestep ;    // Table index
                __m256d tpre, tpim;
                // in GS, calc
                //  x = x + y
                //  y = w(x - y)
                tpre =  _mm256_fmadd_pd (realv[vi + halfsize], trigtablesv[ti], imagv[vi + halfsize] * trigtablesv[ti + 1]);
                tpim =  _mm256_fmsub_pd (imagv[vi + halfsize], trigtablesv[ti], realv[vi + halfsize] * trigtablesv[ti + 1]);
                realv[vi + halfsize] = realv[vi] - tpre;
                imagv[vi + halfsize] = imagv[vi] - tpim;

                realv[vi] += tpre;
                imagv[vi] += tpim;
            }
        }
        if (size == n/4)   // Prevent overflow in 'size *= 2'
            break;

        trigtablesv += size;
    }
#else
    fp_t re, im;

    for (size = 8; size <= n; size <<= 1)
    {
        halfsize = size >> 1;

#ifdef VERBOSE
        printf("\nCT Round #%d: {\n", floor_log2(size));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        for (i = 0; i < n; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 4, tablestep += 8)
            {
                for (size_t k = 0; k < 4; k++)   // To simulate x86 AVX 4-vectors
                {
                    uint64_t vi = i + j + k;        // Vector index
                    uint64_t ti = tablestep + k;    // Table index

                    // in CT, calc
                    //  x = x + wy
                    //  y = x - wy

                    re = real[vi + halfsize];   // y
                    im = imag[vi + halfsize];
                    tpre = re * trigtables[ti] + im * trigtables[ti + 4];   // wy
                    tpim = im * trigtables[ti] - re * trigtables[ti + 4];

                    real[vi + halfsize] = real[vi] - tpre;   // y = x - wy
                    imag[vi + halfsize] = imag[vi] - tpim;
                    real[vi] += tpre;   // x = x + wy
                    imag[vi] += tpim;
                }
            }
        }
        if (size == n)   // Prevent overflow in 'size *= 2'
            break;

        trigtables += size;
    }

#endif

#ifdef VERBOSE
    printf("\nCT Result:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

}


// ---   Implementation of FFT via the Gentleman-Sande data path   -------------

static void dft_GS(const FftTables *const tables,
                   fp_t *const real,
                   fp_t *const imag,
                   const bool inverse,
                   const bool bitrev)
{
    // check table type
    if (!((inverse && tables->type == IFFT_e) || (!inverse && tables->type == FFT_e))) return;

    const size_t n = tables->n;

    size_t i, j;
    fp_t re, im;

    size_t halfsize, tablestep, size;

    // Gentleman-Sande decimation-in-frequency
    const fp_t * trigtables = tables->gs_path_tables;   // no *const since the pointer is iterated

#ifdef FMA_OPT
    __m256d * trigtablesv = (__m256d *) trigtables;
    for (size = n/4; size >= 2; size >>= 1)
    {
        halfsize = size >> 1;

        for (i = 0; i < n/4; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 1, tablestep += 2){
                uint64_t vi = i + j;        // Vector index
                uint64_t ti = tablestep ;    // Table index
                __m256d * realv = (__m256d *) real;
                __m256d * imagv = (__m256d *) imag;
                __m256d rev, imv;
                // in GS, calc
                //  x = x + y
                //  y = w(x - y)
                rev = realv[vi] - realv[vi + halfsize]; // x - y
                imv = imagv[vi] - imagv[vi + halfsize];

                realv[vi] += realv[vi + halfsize];   // x = x + y
                imagv[vi] += imagv[vi + halfsize];

                realv[vi + halfsize] = _mm256_fmadd_pd (rev, trigtablesv[ti], imv * trigtablesv[ti + 1]);
                imagv[vi + halfsize] = _mm256_fmsub_pd (imv, trigtablesv[ti], rev * trigtablesv[ti + 1]);
                // real[vi + halfsize] = re * trigtables[ti] + im * trigtables[ti + 4];   // y = w(x - y)
                // imag[vi + halfsize] = im * trigtables[ti] - re * trigtables[ti + 4];
            }
        }

        trigtablesv += size;
    }
#else
    for (size = n; size >= 8; size >>= 1)
    {
        halfsize = size >> 1;

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(size));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        for (i = 0; i < n; i += size)
        {
            for (j = 0, tablestep = 0;  j < halfsize;   j += 4, tablestep += 8){
                for (size_t k = 0; k < 4; k++)   // To simulate x86 AVX 4-vectors
                {
                    uint64_t vi = i + j + k;        // Vector index
                    uint64_t ti = tablestep + k;    // Table index

                    // in GS, calc
                    //  x = x + y
                    //  y = w(x - y)
                    re = real[vi] - real[vi + halfsize];   // x - y
                    im = imag[vi] - imag[vi + halfsize];

                    real[vi] += real[vi + halfsize];   // x = x + y
                    imag[vi] += imag[vi + halfsize];
                    real[vi + halfsize] = re * trigtables[ti] + im * trigtables[ti + 4];   // y = w(x - y)
                    imag[vi + halfsize] = im * trigtables[ti] - re * trigtables[ti + 4];
                }
            }
        }

        trigtables += size;
    }
#endif

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(4));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 4 merge (special)
    if (n >= 4)
    {
        if (!inverse)
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices (same as CT)
                re = real[i + 2];   // y
                im = imag[i + 2];
                real[i + 2] = real[i] - re;   // y = x - y
                imag[i + 2] = imag[i] - im;
                real[i] += re;   // x = x + y
                imag[i] += im;

                // odd indices (times i)
                re =  (imag[i + 1] - imag[i + 3]);   // -i(x - y)
                im = -(real[i + 1] - real[i + 3]);
                real[i + 1] += real[i + 3];   // x = x + y
                imag[i + 1] += imag[i + 3];
                real[i + 3] = re;   // y = i(x - y)
                imag[i + 3] = im;
            }
        }
        else
        {
            for (i = 0; i < n; i += 4)
            {
                // even indices (same as CT)
                re = real[i + 2];   // y
                im = imag[i + 2];
                real[i + 2] = real[i] - re;   // y = x - y
                imag[i + 2] = imag[i] - im;
                real[i] += re;   // x = x + y
                imag[i] += im;

                // odd indices (times -i)
                re = -(imag[i + 1] - imag[i + 3]);   // i(x - y)
                im =  (real[i + 1] - real[i + 3]);
                real[i + 1] += real[i + 3];   // x = x + y
                imag[i + 1] += imag[i + 3];
                real[i + 3] = re;   // y = i(x - y)
                imag[i + 3] = im;
            }
        }
    }

#ifdef VERBOSE
        printf("\nGS Round #%d: {\n", floor_log2(2));
        for (int ii = 0; ii < n/2; ii++)
            printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

    // size 2 merge (special; same as CT)
    if (n >= 2)
    {
        for (i = 0; i < n; i += 2)
        {
            // addition and subtraction
            re = real[i];
            im = imag[i];

            real[i] += real[i + 1];   // x = x + y
            imag[i] += imag[i + 1];
            real[i + 1] = re - real[i + 1];   // y = x - y
            imag[i + 1] = im - imag[i + 1];
        }
    }

    if (bitrev)
    {
#ifdef VERBOSE
    printf("\nBefore bit-rev:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

        const size_t *const bitreversed = tables->bit_reversed;
        permute_ary(bitreversed, real,
                    imag,
                    n, inverse);
    }
    else if (inverse)   // n.b., division is a part of bitreverse
    {
        for (i = 0; i < n; i++)
        {
            real[i] /= n;
            imag[i] /= n;
        }
    }

#ifdef VERBOSE
    printf("\nGS Result:\n");
    for (int ii = 0; ii < n/2; ii++)
        printf(" %+.3f%+.3fi | %+.3f%+.3fi\n", real[2*ii], imag[2*ii], real[2*ii+1], imag[2*ii+1]);
#endif   // #ifdef VERBOSE

}


// ---   Helper functions   ----------------------------------------------------

static inline void permute_ary(const size_t *const permutation,
                               fp_t *const re,
                               fp_t *const im,
                               const size_t n,
                               const bool inverse)
{
    if (inverse)
    {
        for (size_t i = 0; i < n; i++)
        {
            size_t j = permutation[i];
            if (i <= j)
            {
                fp_t tp0re = re[i];
                fp_t tp1re = re[j];
                re[i] = tp1re / n;
                re[j] = tp0re / n;

                fp_t tp0im = im[i];
                fp_t tp1im = im[j];
                im[i] = tp1im / n;
                im[j] = tp0im / n;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < n; i++)
        {
            size_t j = permutation[i];
            if (i < j)
            {
                fp_t tp0re = re[i];
                fp_t tp1re = re[j];
                re[i] = tp1re;
                re[j] = tp0re;

                fp_t tp0im = im[i];
                fp_t tp1im = im[j];
                im[i] = tp1im;
                im[j] = tp0im;
            }
        }
    }
}

// n must be a multiple of 4
static fp_t accurate_sine(uint64_t i, uint64_t n)
{
    if (n % 4 != 0)
    {
        return (fp_t)NAN;
    }
    else
    {
        int32_t neg = 0;  // Boolean
        // Reduce to full cycle
        i %= n;
        // Reduce to half cycle
        if (i >= n / 2)
        {
            neg = 1;
            i -= n / 2;
        }
        // Reduce to quarter cycle
        if (i >= n / 4)
            i = n / 2 - i;
        // Reduce to eighth cycle
        fp_t val;
        if (i * 8 < n)
            val = (fp_t)sinl(2 * (fp_t)M_PI * i / n);
        else
            val = (fp_t)cosl(2 * (fp_t)M_PI * (n / 4 - i) / n);
        // Apply sign
        return neg ? -val : val;
    }
}

static int32_t floor_log2(size_t n)
{
    int32_t result = 0;
    for (; n > 1; n /= 2)
        result++;
    return result;
}

static size_t reverse_bits(size_t x, uint32_t n)
{
    size_t result = 0;
    uint32_t i;
    for (i = 0; i < n; i++, x >>= 1)
        result = (result << 1) | (x & 1);
    return result;
}


// ---   Destroy tables   ------------------------------------------------------

void tables_destroy(FftTables *const tables)
{
    if (tables == NULL)
        return;
    free(tables->bit_reversed);
    free(tables->cos_table);
    free(tables->sin_table);
    free(tables->ct_path_tables);
    free(tables->gs_path_tables);
    memset(tables, 0, sizeof(FftTables));
    free(tables);
}


/* Functions for TFHE */

FFT_Processor_FFNT new_FFT_Processor_FFNT(int N){
    FFT_Processor_FFNT res;
    res = (FFT_Processor_FFNT) safe_malloc(sizeof(*res));
    res->N = N;
    res->fpr = (double *) safe_aligned_malloc(sizeof(double) * N);
    res->fpi = (double *) safe_aligned_malloc(sizeof(double) * N);
    res->fft_n_2_tables = fft_init(N/2);
    res->ifft_n_2_tables = ifft_init(N/2);
    res->ffnt_2n_tables  = fft_init(2*N);
    return res;
}

void execute_reverse_torus64(double * res, const uint64_t * a, FFT_Processor_FFNT proc){
    const int N = proc->N;
    int64_t * aa = (int64_t *) a;
    for (size_t i = 0; i < N; i++){
        res[i] = ((double) aa[i]);
    }
    ffnt_transform(proc->ffnt_2n_tables, proc->fft_n_2_tables, res, proc->fpi);
    memcpy(&res[N/2], proc->fpi, sizeof(double) * N/2);
    // for (size_t i = 0; i < N/2; i++){
    //     res[i + N/2] = proc->fpi[i];
    // }
}

void execute_direct_torus64(uint64_t * res, const double * a, FFT_Processor_FFNT proc){
    const int N = proc->N;
    for (size_t i = 0; i < N/2; i++){
        proc->fpr[i] = a[i]; 
        proc->fpi[i] = a[i + N/2];
    }
    memset(&proc->fpr[N/2], 0, sizeof(double)*N/2);
    memset(&proc->fpi[N/2], 0, sizeof(double)*N/2);
    // for (size_t i = N/2; i < N; i++){
    //     proc->fpr[i] = 0.; proc->fpi[i] = 0.;
    // }
    iffnt_transform(proc->ffnt_2n_tables, proc->ifft_n_2_tables, proc->fpr, proc->fpi);
    const uint64_t* const vals = (const uint64_t*) proc->fpr;
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    for (size_t i = 0; i < N; i++){
        uint64_t val = (vals[i]&valmask0)|valmask1; //mantissa on 53 bits
        uint16_t expo = (vals[i]>>52)&expmask0; //exponent 11 bits
        // 1023 -> 52th pos -> 0th pos
        // 1075 -> 52th pos -> 52th pos
        int16_t trans = expo-1075;
        uint64_t val2 = trans>0?(val<<trans):(val>>-trans);
        res[i] = (vals[i]>>63)?-val2:val2;
    }
}


#ifdef __cplusplus
}
#endif
