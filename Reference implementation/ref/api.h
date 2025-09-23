#ifndef API_H
#define API_H

#include "params.h"

#define CRYPTO_SECRETKEYBYTES  LORE_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES  LORE_PUBLICKEYBYTES
#define CRYPTO_CIPHERTEXTBYTES LORE_CIPHERTEXTBYTES
#define CRYPTO_BYTES           LORE_SYMBYTES

#if   (LORE_K == 2)
#define CRYPTO_ALGNAME "Lore-128"
#elif (LORE_K == 3)
#define CRYPTO_ALGNAME "Lore-192"
#elif (LORE_K == 4)
#define CRYPTO_ALGNAME "Lore-256"
#endif

// PKE API
int crypto_pke_keypair(unsigned char *pk, unsigned char *sk);

int crypto_pke_encrypt(unsigned char *ct,
                       const unsigned char *m,
                       unsigned long long mlen,
                       const unsigned char *pk);

int crypto_pke_decrypt(unsigned char *m,
                       unsigned long long *mlen,
                       const unsigned char *ct,
                       const unsigned char *sk);

#endif