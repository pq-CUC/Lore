#ifndef PKE_H
#define PKE_H

#include "params.h"

// These are wrappers for the API functions, so they don't need namespace prefixes.
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