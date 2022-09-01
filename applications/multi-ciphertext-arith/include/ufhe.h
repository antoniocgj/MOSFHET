#include <mosfhet.h>
#include <stdbool.h>

typedef struct{
  TLWE_Key tlwe;
  TRLWE_Key trlwe;
  TLWE_Key extracted_trlwe;
  TRGSW_Key trgsw;
} * ufhe_priv_keyset;


typedef struct{
  Bootstrap_Key bootstrap_key;
  LUT_Packing_KS_Key packing_key;
  TLWE_KS_Key ks_key;
} * ufhe_public_keyset;

typedef struct{
  TLWE * digits;
  int d;
  bool _signed;
} * ufhe_integer;

struct aux_data{
  Torus * lut_torus_slots;
  TLWE * lut_tlwe_slots;
  TLWE tlwe_zero, tlwe_selector, * tlwe_lut;
  TRLWE ADDSUB_LUT, SIGNEXTEND_LUT, COMPARE_LUT, * LUT;
  int ** mulmod_matrix, ** mulquo_matrix;
};


typedef struct{
  int torus_base, log_torus_base;
  ufhe_public_keyset keyset;
  struct aux_data aux;

} * ufhe_context;

// TODO: 
typedef enum{
  SET0,
  SET1,
  SET2
} ufhe_params;


/* Setup */
ufhe_priv_keyset ufhe_new_priv_keyset(ufhe_params params);
ufhe_public_keyset ufhe_new_public_keyset(ufhe_priv_keyset priv_key, ufhe_params params);
ufhe_context ufhe_setup_context(ufhe_public_keyset keyset);

/* IO */
ufhe_priv_keyset ufhe_load_priv_keyset(FILE * fd);
ufhe_public_keyset ufhe_load_public_keyset(FILE * fd);
void ufhe_save_priv_keyset(FILE * fd, ufhe_priv_keyset keyset);
void ufhe_save_public_keyset(FILE * fd, ufhe_public_keyset keyset);

/* Integer */
ufhe_integer ufhe_new_integer(int precision, bool _signed, ufhe_context ctx);
void ufhe_free_integer(ufhe_integer c);
void ufhe_cleartext_integer(ufhe_integer out, uint64_t value, ufhe_context ctx);
void ufhe_encrypt_integer(ufhe_integer out, uint64_t value, ufhe_priv_keyset key, ufhe_context ctx);
int64_t ufhe_decrypt_integer(ufhe_integer c, ufhe_priv_keyset key, ufhe_context ctx);

/* Arithmetic */

void ufhe_extend_integer(ufhe_integer c, int old_precision, ufhe_context ctx);

/* b = a */
void ufhe_copy_integer(ufhe_integer b, ufhe_integer a, ufhe_context ctx);

/* c = a + b */
void ufhe_add_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx);

/* c = a*B^g + b*B^h */
void ufhe_sl_add_integer(ufhe_integer c, ufhe_integer a, int g, ufhe_integer b, int h, ufhe_context ctx);

/* c = a - b */
void ufhe_sub_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx);

/* c = a * b */
void ufhe_mul_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx);

/* Logic */

/* c =     
        0, if a < b, 
        1, if a = b,
        2, if a > b */
void ufhe_cmp_integer(ufhe_integer c, ufhe_integer a, ufhe_integer b, ufhe_context ctx);

/* out = {a, b, ...}[selector] */
void ufhe_mux_integer(ufhe_integer out, ufhe_integer selector, ufhe_context ctx, ufhe_integer a, ufhe_integer b, ...);

/* out = vec[selector] */
void ufhe_mux_integer_array(ufhe_integer out, ufhe_integer selector, ufhe_context ctx, int size, ufhe_integer * vec);

/* out = lut[selector] */
void ufhe_lut_integer(ufhe_integer out, ufhe_integer selector, uint64_t * lut, int size, ufhe_context ctx);

/* ML functions */

/* out = in > 0 ? in : 0 */
void ufhe_relu_integer(ufhe_integer out, ufhe_integer in, ufhe_context ctx);