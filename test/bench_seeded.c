#include "mosfhet.h"
#include <sys/time.h>
#include <time.h>

uint64_t get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_usec) + (tv.tv_sec * 1000000);
}
static uint64_t ts_ns()
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000LL + (uint64_t)ts.tv_nsec;
}


uint64_t __g_clock_begin, __g_clock_end;
#define MEASURE_TIME(REP, MSG, CODE) \
  __g_clock_begin = ts_ns(); \
  for (size_t __i = 0; __i < REP; __i++){\
    CODE;\
  }\
  __g_clock_end = ts_ns(); \
  printf("%lu%03lu%03lu", ((__g_clock_end - __g_clock_begin)/REP)/1000000, (((__g_clock_end - __g_clock_begin)/REP)/1000)%1000, ((__g_clock_end - __g_clock_begin)/REP)%1000);
  // printf(MSG ": %lu,%03lu,%03lu ns\n", ((__g_clock_end - __g_clock_begin)/REP)/1000000, (((__g_clock_end - __g_clock_begin)/REP)/1000)%1000, ((__g_clock_end - __g_clock_begin)/REP)%1000);

// LWE params
const int n = 630;
const double lwe_std_dev = 3.0517578125e-05; // 2^-15
// RLWE params
const int N = 2048, k = 1;
const double rlwe_std_dev = 5.684341886080802e-14; // 2^-44
// RGSW params
const int l = 6, Bg_bit = 7;
// KS params
const int t = 6, base_bit = 2;

#define PRINT_SIZE 20

void print_poly_uint64(TorusPolynomial p){
  for (size_t i = 0; i < PRINT_SIZE; i++){
    printf("%lx, ", p->coeffs[i]);
  }
  printf("\n");
}

void print_poly_double(DFT_Polynomial p){
  for (size_t i = 0; i < PRINT_SIZE; i++){
    printf("%lf, ", p->coeffs[i]);
  }
  printf("\n");
}

void print_poly_double_format(DFT_Polynomial p){
  const uint64_t exp_mask = 0x7FF;
  uint64_t * c = (uint64_t *) p->coeffs;
  for (size_t i = 0; i < PRINT_SIZE; i++){
    printf("%ld, ", (c[i]>>52)&exp_mask);
  }
  printf("\n");
  for (size_t i = 0; i < PRINT_SIZE; i++){
    printf("%ld, ", (((c[i]>>52)&exp_mask) - 1023));
  }
  printf("\n");
}

#define SAMPLES 5000
#define __REP (SAMPLES*100)

int main(int argc, char const *argv[]){
  TRLWE_Key key = trlwe_new_binary_key(N, k, rlwe_std_dev);
  TRLWE out = trlwe_new_sample(0, key);
  // DFT_Polynomial p = polynomial_new_DFT_polynomial(N);
  // polynomial_torus_to_DFT(p, out->a[0]);
  // print_poly_double_format(p);

  TRLWE c[SAMPLES], c_cmp[SAMPLES];
  for (size_t i = 0; i < SAMPLES; i++){
    c_cmp[i] = trlwe_new_compressed_sample(0, key);
    c[i] = trlwe_new_sample(0, key);
  }
  
  // init_fft(N);

  for (size_t i = 0; i < 4000; i+=50){
    uint64_t j = i; if(j == 0) j = 1;

    MEASURE_TIME(__REP, "TRLWE Subto ", 
      trlwe_subto(out, c[__i%j]);
    );
    printf("\t");
    MEASURE_TIME(__REP, "TRLWE comp. subto: ", 
      trlwe_compressed_subto(out, c_cmp[__i%j]);
    );
    printf("\n");
  }
  

  // MEASURE_TIME(__REP, "TRLWE Subto ", 
  //   trlwe_subto(out, c[__i%SAMPLES]);
  // );

  // MEASURE_TIME(__REP, "TRLWE comp. subto: ", 
  //   trlwe_compressed_subto(out, c_cmp[__i%SAMPLES]);
  // );

  return 0;
}
