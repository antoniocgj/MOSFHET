#include "mosfhet.h"

double torus2double(Torus x){
  return ((double) x)/18446744073709551615.0;
}

Torus double2torus(double x){
  return (Torus) ((int64_t) (((double)18446744073709551615.0)*x));
}

/* return round(in[i] * 2^log_scale) */
uint64_t torus2int(Torus x, int log_scale){
  const uint64_t bit_size = sizeof(Torus) * 8;
  const Torus round_offset = 1UL << (bit_size - log_scale - 1);
  return (x + round_offset)>>(bit_size - log_scale);
}

/* return round(in[i] / 2^log_scale) */
Torus int2torus(uint64_t x, int log_scale){
  const uint64_t bit_size = sizeof(Torus) * 8;
  return x << (bit_size - log_scale);
}

// Random generation

#ifndef PORTABLE_BUILD
// TODO: add code src.
void generate_rnd_seed(uint64_t * p){
  if(0 == _rdrand64_step ((unsigned long long *) p) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[1])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[2])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[3]))){
    printf("Random Generation Failed\n");
    return;
  }
}
#else 
void generate_rnd_seed(uint64_t * p){
  FILE *fp;
  fp = fopen("/dev/urandom", "r");
  fread(p, 1, 32, fp);
  fclose(fp);
}
#endif


#include "sha3/fips202.h"

void get_rnd_from_hash(uint64_t amount, uint8_t * pointer){
  uint64_t rnd[4];
  generate_rnd_seed(rnd);
  shake256(pointer, amount, (uint8_t *) rnd, 32);
}

void get_rnd_from_buffer(uint64_t amount, uint8_t * pointer){
  static uint8_t buffer[1024];
  static int idx = 1024;
  if(amount > (1024 - idx)){
    idx = 0;
    get_rnd_from_hash(1024, buffer);
  }
  memcpy(pointer, buffer+idx, amount);
  idx += amount;
}

void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  if(amount < 512) get_rnd_from_buffer(amount, pointer);
  else get_rnd_from_hash(amount, pointer);
}

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
double generate_normal_random(double sigma){
  Torus rnd[2];
  generate_random_bytes(16, (uint8_t *) rnd);
  return cos(2.*M_PI*torus2double(rnd[0]))*sqrt(-2.*log(torus2double(rnd[1])))*sigma;
}

void generate_torus_normal_random_array(Torus * out, double sigma, int N){
  for (size_t i = 0; i < N; i++){
    out[i] = double2torus(generate_normal_random(sigma));
  }
}

// Mem alloc

uint64_t _glb_mem_count = 0;
// safe_malloc
// https://stackoverflow.com/questions/48043811/creating-a-function-to-check-if-malloc-succeeded
void * safe_malloc(size_t size){
  void *ptr = malloc(size);
  if (!ptr && (size > 0)) {
    perror("malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size); 
  return ptr;
}

void * safe_aligned_malloc(size_t size){
  void * ptr;
  #ifdef AVX512_OPT
  int err = posix_memalign(&ptr, 64, size);
  #else
  int err = posix_memalign(&ptr, 32, size);
  #endif
  if (err || (!ptr && (size > 0))) {
    perror("aligned malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size);
  return ptr;
}