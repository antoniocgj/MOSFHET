#include "unity_test/unity.h"
#include "ufhe.h"

void setUp(void) {}
void tearDown(void) {}


ufhe_priv_keyset _priv_key;
ufhe_public_keyset _pub_key;
ufhe_context _ctx;

void test_setup(){
  _priv_key = ufhe_new_priv_keyset(SET0);
  _pub_key = ufhe_new_public_keyset(_priv_key, SET0);
  _ctx = ufhe_setup_context(_pub_key);
}

void test_setup_save_and_exit(){
  test_setup();
  FILE * pub_fd = fopen("pub.key", "w");
  ufhe_save_public_keyset(pub_fd, _pub_key);
  fclose(pub_fd);
  FILE * priv_fd = fopen("priv.key", "w");
  ufhe_save_priv_keyset(priv_fd, _priv_key);
  fclose(priv_fd);
  exit(0);
}

void test_setup_from_io(){
  FILE * priv_fd = fopen("priv.key", "r");
  _priv_key = ufhe_load_priv_keyset(priv_fd);
  fclose(priv_fd);
  FILE * pub_fd = fopen("pub.key", "r");
  _pub_key = ufhe_load_public_keyset(pub_fd);
  fclose(pub_fd);
  _ctx = ufhe_setup_context(_pub_key);
}

void test_int_encrypt(){
  ufhe_integer c = ufhe_new_integer(8, true, _ctx);
  int8_t ct;
  for (size_t i = 0; i < 1000; i++){
    generate_random_bytes(1, (uint8_t *) &ct);
    ufhe_encrypt_integer(c, ct, _priv_key, _ctx);
    int8_t res = ufhe_decrypt_integer(c, _priv_key, _ctx);
    TEST_ASSERT_EQUAL_INT8_MESSAGE(ct, res, "Integer encryption failed");
  }
  
}

void test_int_extend(){
  ufhe_integer c = ufhe_new_integer(32, true, _ctx);
  int8_t ct;
  generate_random_bytes(1, (uint8_t *) &ct);
  ct |= 0x80;

  c->d = 8/_ctx->log_torus_base;
  ufhe_encrypt_integer(c, ct, _priv_key, _ctx);
  c->d = 32/_ctx->log_torus_base;
  ufhe_extend_integer(c, 8, _ctx);

  int32_t res = ufhe_decrypt_integer(c, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT32_MESSAGE((int32_t) ct, res, "Integer extend negative failed");
  
  generate_random_bytes(1, (uint8_t *) &ct);
  ct &= 0x80 - 1;

  c->d = 8/_ctx->log_torus_base;
  ufhe_encrypt_integer(c, ct, _priv_key, _ctx);
  c->d = 32/_ctx->log_torus_base;
  ufhe_extend_integer(c, 8, _ctx);

  res = ufhe_decrypt_integer(c, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT32_MESSAGE((int32_t) ct, res, "Integer extend positive failed");
}

void test_int_add(){
  ufhe_integer a = ufhe_new_integer(8, true, _ctx);
  ufhe_integer b = ufhe_new_integer(8, true, _ctx);
  ufhe_integer c = ufhe_new_integer(8, true, _ctx);
  int8_t ct[2];
  for (size_t i = 0; i < 100; i++){
    generate_random_bytes(2, (uint8_t *) ct);
    ufhe_encrypt_integer(a, ct[0], _priv_key, _ctx);
    ufhe_encrypt_integer(b, ct[1], _priv_key, _ctx);
    ufhe_add_integer(c, a, b, _ctx);
    int8_t res = ufhe_decrypt_integer(c, _priv_key, _ctx);
    TEST_ASSERT_EQUAL_INT8_MESSAGE((ct[0] + ct[1])%256, res, "Integer add failed");

    ufhe_sub_integer(c, a, b, _ctx);
    res = ufhe_decrypt_integer(c, _priv_key, _ctx);
    TEST_ASSERT_EQUAL_INT8_MESSAGE((ct[0] - ct[1])%256, res, "Integer sub failed");
  }  
}


void test_int_sl_add(){
  const int log_B = _ctx->log_torus_base, mask = 0xFF;
  ufhe_integer a = ufhe_new_integer(8, true, _ctx);
  ufhe_integer b = ufhe_new_integer(8, true, _ctx);
  ufhe_integer c = ufhe_new_integer(32, true, _ctx);
  int8_t ct[2];
  generate_random_bytes(2, (uint8_t *) ct);
  ct[0] = 125;
  ct[1] = -22;
  for (size_t i = 0; i < 4; i++){
    for (size_t j = 0; j < 4; j++){
      ufhe_encrypt_integer(a, ct[0], _priv_key, _ctx);
      ufhe_encrypt_integer(b, ct[1], _priv_key, _ctx);
      ufhe_sl_add_integer(c, a, i, b, j, _ctx);
      int res = ufhe_decrypt_integer(c, _priv_key, _ctx);
      TEST_ASSERT_EQUAL_INT8_MESSAGE(((ct[0]<<(log_B*i)) + (ct[1]<<(log_B*j)))&mask, res, "Integer sl add failed");
    }
  }
}

void test_int_mul(){
  const int mask = 0xFF;
  ufhe_integer a = ufhe_new_integer(8, true, _ctx);
  ufhe_integer b = ufhe_new_integer(8, true, _ctx);
  ufhe_integer c = ufhe_new_integer(32, true, _ctx);
  int8_t ct[2];
  generate_random_bytes(2, (uint8_t *) ct);
  ufhe_encrypt_integer(a, ct[0], _priv_key, _ctx);
  ufhe_encrypt_integer(b, ct[1], _priv_key, _ctx);
  ufhe_mul_integer(c, a, b, _ctx);
  int res = ufhe_decrypt_integer(c, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT8_MESSAGE((ct[0]*ct[1])&mask, res, "Signed Integer mul failed");

  a->_signed = false;
  b->_signed = false;
  c->_signed = false;

  uint8_t ct2[2];
  generate_random_bytes(2, (uint8_t *) ct2);
  ufhe_encrypt_integer(a, ct2[0], _priv_key, _ctx);
  ufhe_encrypt_integer(b, ct2[1], _priv_key, _ctx);
  ufhe_mul_integer(c, a, b, _ctx);
  uint32_t res2 = ufhe_decrypt_integer(c, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT32_MESSAGE(ct2[0]*ct2[1], res2, "Unsigned Integer mul failed");
}

void test_int_cmp(){
  ufhe_integer a = ufhe_new_integer(8, true, _ctx);
  ufhe_integer b = ufhe_new_integer(8, true, _ctx);
  ufhe_integer c = ufhe_new_integer(8, true, _ctx);
  int8_t ct[2];
  uint8_t ct2[2];
  uint32_t res2;
  int res, expc;

  for (size_t i = 0; i < 100; i++){
    a->_signed = true;
    b->_signed = true;
    c->_signed = true;

    generate_random_bytes(2, (uint8_t *) ct);
    ufhe_encrypt_integer(a, ct[0], _priv_key, _ctx);
    ufhe_encrypt_integer(b, ct[1], _priv_key, _ctx);
    ufhe_cmp_integer(c, a, b, _ctx);
    res = ufhe_decrypt_integer(c, _priv_key, _ctx);
    expc = ct[0]>ct[1]? 2 : ct[0] == ct[1] ? 1 : 0;
    TEST_ASSERT_EQUAL_INT8_MESSAGE(expc, res, "Signed Integer cmp failed");

    a->_signed = false;
    b->_signed = false;
    c->_signed = false;

    generate_random_bytes(2, (uint8_t *) ct2);
    ufhe_encrypt_integer(a, ct2[0], _priv_key, _ctx);
    ufhe_encrypt_integer(b, ct2[1], _priv_key, _ctx);
    ufhe_cmp_integer(c, a, b, _ctx);
    res2 = ufhe_decrypt_integer(c, _priv_key, _ctx);
    expc = ct2[0]>ct2[1]? 2 : ct2[0] == ct2[1] ? 1 : 0;
    TEST_ASSERT_EQUAL_INT32_MESSAGE(expc, res2, "Unsigned Integer cmp failed");
  }
  
}

void test_relu(){
  ufhe_integer a = ufhe_new_integer(8, true, _ctx);
  int8_t ct, res;
  for (size_t i = 0; i < 10; i++){
    generate_random_bytes(1, (uint8_t *) &ct);
    ufhe_encrypt_integer(a, ct, _priv_key, _ctx);
    ufhe_relu_integer(a, a, _ctx);
    res = ufhe_decrypt_integer(a, _priv_key, _ctx);
    TEST_ASSERT_EQUAL_INT8_MESSAGE((ct > 0? ct : 0), res, "ReLU failed");
  }
}



void test_lut(){
  ufhe_integer selector = ufhe_new_integer(4, true, _ctx);
  ufhe_integer out = ufhe_new_integer(8, true, _ctx);
  int8_t ct, res, ct_lut[16];
  ufhe_integer vec[16];
  generate_random_bytes(16, (uint8_t *) ct_lut);
  generate_random_bytes(1, (uint8_t *) &ct);
  ct &= 0xf;
  ufhe_encrypt_integer(selector, ct, _priv_key, _ctx);
  for (size_t i = 0; i < 16; i++){
    vec[i] = ufhe_new_integer(8, true, _ctx);
    ufhe_encrypt_integer(vec[i], ct_lut[i], _priv_key, _ctx);
  }
  ufhe_mux_integer_array(out, selector, _ctx, 16, vec);
  res = ufhe_decrypt_integer(out, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT8_MESSAGE(ct_lut[ct], res, "MUX failed");
}


void test_lut_ct(){
  ufhe_integer selector = ufhe_new_integer(4, true, _ctx);
  ufhe_integer out = ufhe_new_integer(8, true, _ctx);
  int8_t ct, res;
  uint64_t ct_lut[16];
  generate_random_bytes(16*8, (uint8_t *) ct_lut);
  generate_random_bytes(1, (uint8_t *) &ct);
  ct &= 0xf;
  for (size_t i = 0; i < 16; i++) ct_lut[i] &= 0xff;
  ufhe_encrypt_integer(selector, ct, _priv_key, _ctx);
  ufhe_lut_integer(out, selector, ct_lut, 16, _ctx);
  res = ufhe_decrypt_integer(out, _priv_key, _ctx);
  TEST_ASSERT_EQUAL_INT8_MESSAGE(ct_lut[ct], res, "LUT failed");
}


int main(int argc, char const *argv[])
{
  // test_setup_save_and_exit();
  test_setup();
  // test_setup_from_io();
  UNITY_BEGIN();
  RUN_TEST(test_lut_ct);
  RUN_TEST(test_lut);
  RUN_TEST(test_relu);
  RUN_TEST(test_int_encrypt);
  RUN_TEST(test_int_extend);
  RUN_TEST(test_int_sl_add);
  RUN_TEST(test_int_cmp);
  RUN_TEST(test_int_mul);
  RUN_TEST(test_int_add);
  return UNITY_END();
}