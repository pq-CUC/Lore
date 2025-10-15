#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../api.h"
#include "../indcpa.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../randombytes.h"
#include "../randombytes.h"
#include "../symmetric.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 1000
#define MLEN CRYPTO_BYTES

uint64_t t[NTESTS];

// Static variables to avoid stack overflow
static poly_crt_vec a_row;
static poly_crt_vec b_dense;
static poly_crt_vec b_crt_sparse;
static poly_sparse b_t_sparse[LORE_K];
static poly_crt res_crt;
// ==================================

/*************************************************
* Name:        main
*
* Description: Main function for the speed test program.
* It benchmarks the performance of the cryptographic operations.
*
* Arguments:   - void
*
* Returns 0 (success)
**************************************************/
int main(void)
{
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sk_pke[CRYPTO_PKE_SECRETKEYBYTES];
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
  unsigned char msg[MLEN];
  unsigned long long mlen;
  int i;

  unsigned char sk_kem[CRYPTO_SECRETKEYBYTES];
  unsigned char ss[CRYPTO_BYTES];
  unsigned char coins_kp[CRYPTO_KEYPAIRCOINBYTES];
  randombytes(coins_kp, CRYPTO_KEYPAIRCOINBYTES);

  
  // Variables for micro-benchmarking
  unsigned char seed[LORE_SYMBYTES];
  unsigned char coins[LORE_SYMBYTES];
  poly_crt_vec a[LORE_K];
  poly_crt_vec sp_crt;
  poly_sparse sp_sparse[LORE_K];


  randombytes(seed, LORE_SYMBYTES);
  randombytes(coins, LORE_SYMBYTES);

  printf("\nBenchmarking %s\n\n", CRYPTO_ALGNAME);

  // === Top Level API ===
  printf("--- Top Level API ---\n");
  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      crypto_pke_keypair(pk, sk_pke);
  }
  print_results("crypto_pke_keypair:", t, NTESTS);

  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      crypto_pke_encrypt(ct, msg, MLEN, pk);
  }
  print_results("crypto_pke_encrypt:", t, NTESTS);

  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      crypto_pke_decrypt(msg, &mlen, ct, sk_pke);
  }
  print_results("crypto_pke_decrypt:", t, NTESTS);


  // === IND-CPA Layer ===
  printf("\n--- IND-CPA Layer ---\n");
  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      indcpa_keypair(pk, sk_pke);
  }
  print_results("indcpa_keypair:", t, NTESTS);

  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      indcpa_enc(ct, msg, pk, coins);
  }
  print_results("indcpa_enc:", t, NTESTS);

  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      indcpa_dec(msg, ct, sk_pke);
  }
  print_results("indcpa_dec:", t, NTESTS);


  // === Polynomial Arithmetic & Sampling ===
  printf("\n--- Polynomial Arithmetic & Sampling ---\n");
  // ******** Start of fix ********
  // Test the new, one-step sampling function
  for (i = 0; i < NTESTS; i++) {
    t[i] = cpucycles();
    poly_getnoise(&sp_crt, sp_sparse, coins, (unsigned char)i);
  }
  print_results("poly_getnoise (optimized fixed-weight sampling):", t, NTESTS);

  // NTT/INVNTT tests remain unchanged
  for (i = 0; i < NTESTS; i++) {
    t[i] = cpucycles();
    poly_ntt(&sp_crt.vec[0].q_poly);
  }
  print_results("poly_ntt (q-part):", t, NTESTS);

  for (i = 0; i < NTESTS; i++) {
    t[i] = cpucycles();
    poly_invntt_tomont(&sp_crt.vec[0].q_poly);
  }
  print_results("poly_invntt_tomont (q-part):", t, NTESTS);

  // --- New Benchmarks ---
  // Prepare "dense" vectors a_row and b_dense
  poly_crt_vec temp_matrix[LORE_K];
  gen_matrix_ntt(temp_matrix, seed,0);

  // Take two rows from the generated matrix as our test vectors a_row and b_dense
  // (Assuming LORE_K >= 2)
  a_row = temp_matrix[0];
  b_dense = temp_matrix[1];
  // Prepare "sparse" vectors b_crt_sparse and b_t_sparse
  poly_getnoise(&b_crt_sparse, b_t_sparse, coins, 0);

  // Convert all q-parts to NTT domain
  poly_crt_vec_ntt(&a_row);
  poly_crt_vec_ntt(&b_dense);
  poly_crt_vec_ntt(&b_crt_sparse);

  // Test Karatsuba (DENSE * DENSE)
  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      poly_crt_vec_pointwise_acc_montgomery(&res_crt, &a_row, &b_dense);
  }
  print_results("poly_crt_vec_acc (Karatsuba, DENSE*DENSE):", t, NTESTS);

  // Test our final optimized sparse multiplication (DENSE * SPARSE)
  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      // Use the correct, updated function signature
      poly_crt_vec_pointwise_acc_montgomery_sparse(&res_crt, &a_row, &b_crt_sparse, b_t_sparse);
  }
  print_results("poly_crt_vec_acc_sparse (Optimized, DENSE*SPARSE):", t, NTESTS);
  // ******** End of fix ********


  // === Symmetric Primitives Benchmarking ===
  printf("\n--- Symmetric Primitives ---\n");
  for (i = 0; i < NTESTS; i++) {
      t[i] = cpucycles();
      gen_matrix_std(a, seed);
  }
  print_results("gen_matrix (XOF):", t, NTESTS);

    // === IND-CCA2 KEM Layer ===
    printf("\n--- IND-CCA2 KEM Functions ---\n");
    
    // KEM-specific variables

    randombytes(coins_kp, CRYPTO_KEYPAIRCOINBYTES);

    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_kem_keypair_derand(pk, sk_kem, coins_kp);
    }
    print_results("crypto_kem_keypair_derand:", t, NTESTS);

    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_kem_keypair(pk, sk_kem);
    }
    print_results("crypto_kem_keypair:", t, NTESTS);

    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_kem_enc_derand(ct, ss, pk, coins);
    }
    print_results("crypto_kem_enc_derand:", t, NTESTS);

    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_kem_enc(ct, ss, pk);
    }
    print_results("crypto_kem_enc:", t, NTESTS);

    for (i = 0; i < NTESTS; i++) {
        t[i] = cpucycles();
        crypto_kem_dec(ss, ct, sk_kem);
    }
    print_results("crypto_kem_dec:", t, NTESTS);
  return 0;
}