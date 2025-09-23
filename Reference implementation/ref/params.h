#ifndef PARAMS_H
#define PARAMS_H

// Security level for Lore (1, 2, or 3)
// 1: 128-bit security
// 2: 192-bit security
// 3: 256-bit security
#ifndef LORE_LEVEL
    #error "LORE_LEVEL is not defined. Please specify a level during compilation, e.g., -DLORE_LEVEL=1"
#endif

// Namespace-prefix for all public functions
#define LORE_NAMESPACE(s) pqcrystals_lore_ref_##s

// Common parameters
#define LORE_N 256
#define LORE_Q 257
#define LORE_SYMBYTES 32 // size in bytes of hashes, seeds, etc.

#if LORE_LEVEL == 1
    #define LORE_K 2
    #define LORE_T 3
    #define LORE_R_BITS 0 // |R| = 1 -> log2|R| = 0
    // Hamming weight definitions
    #define LORE_HWT_P1 40
    #define LORE_HWT_M1 40
    #define LORE_HWT_P2 0
    #define LORE_HWT_M2 0
#elif LORE_LEVEL == 2
    #define LORE_K 3
    #define LORE_T 7
    #define LORE_R_BITS 1 // |R| = 2 -> log2|R| = 1
    // Hamming weight definitions
    #define LORE_HWT_P1 50
    #define LORE_HWT_M1 50
    #define LORE_HWT_P2 1
    #define LORE_HWT_M2 1
#elif LORE_LEVEL == 3
    #define LORE_K 4
    #define LORE_T 13
    #define LORE_R_BITS 2 // |R| = 4 -> log2|R| = 2
    // Hamming weight definitions
    #define LORE_HWT_P1 40
    #define LORE_HWT_M1 40
    #define LORE_HWT_P2 22
    #define LORE_HWT_M2 22
#else
    #error "LORE_LEVEL must be defined as 1, 2, or 3"
#endif

// Total Hamming weight
#define LORE_HWT_TOTAL (LORE_HWT_P1 + LORE_HWT_M1 + LORE_HWT_P2 + LORE_HWT_M2)

#if LORE_R_BITS == 0
#define LORE_R_SIZE 1
#elif LORE_R_BITS == 1
#define LORE_R_SIZE 2
#elif LORE_R_BITS == 2
#define LORE_R_SIZE 4
#endif

// Size of the secret key in bytes
#define LORE_SECRETKEYBYTES ( (LORE_K * 2) + (LORE_K * LORE_HWT_TOTAL * 2) + (LORE_K * LORE_POLY_BYTES) )

// Size of the compressed t-part of a polynomial in bytes
#define LORE_POLY_COMPRESSED_BYTES_T ((LORE_N * LORE_R_BITS + 7) / 8)

// Size of the public key in bytes
#define LORE_PUBLICKEYBYTES (LORE_SYMBYTES + (LORE_K * LORE_N) + ((LORE_K * LORE_N + 7) / 8) + (LORE_K * LORE_POLY_COMPRESSED_BYTES_T) + 50)

// Size of the ciphertext in bytes
#define LORE_CIPHERTEXTBYTES (((LORE_K + 1) * LORE_N) + (((LORE_K + 1) * LORE_N + 7) / 8) + (LORE_K * LORE_POLY_COMPRESSED_BYTES_T) + 50)

#define LORE_POLY_BYTES (LORE_N * 2)

#endif