#pragma once

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <mosfhet.h>
extern void * safe_malloc(size_t size);
extern void * safe_aligned_malloc(size_t size);
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

void *new_fft_table(int32_t nn);
double *fft_table_get_buffer(const void *tables);
void *new_ifft_table(int32_t nn);
double *ifft_table_get_buffer(const void *tables);
void fft_model(const void *tables);
void ifft_model(void *tables);
void fft(const void *tables, double *data);
void ifft(const void *tables, double *data);


typedef struct{
    int32_t _2N;
    int32_t N;
    int32_t Ns2;
    double *real_inout_direct;
    double *imag_inout_direct;
    double *real_inout_rev;
    double *imag_inout_rev;
    void *tables_direct;
    void *tables_reverse;
    double *cosomegaxminus1;
    double *sinomegaxminus1;
    int32_t *reva;
} * FFT_Processor_Spqlios;

FFT_Processor_Spqlios new_FFT_Processor_Spqlios(const int32_t N);

void execute_reverse_int(double *res, const int32_t *a, FFT_Processor_Spqlios proc);
void execute_reverse_torus32(double *res, const uint32_t *a, FFT_Processor_Spqlios proc);
void execute_direct_torus32(uint32_t *res, const double *a, FFT_Processor_Spqlios proc);
void execute_reverse_torus64(double* res, const uint64_t* a, FFT_Processor_Spqlios proc);
void execute_direct_torus64(uint64_t* res, const double* a, FFT_Processor_Spqlios proc);
