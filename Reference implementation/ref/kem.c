#include <stddef.h>
#include <string.h>
#include "kem.h"
#include "indcpa.h"
#include "symmetric.h"
#include "randombytes.h"
#include "verify.h"
#include "params.h"
#include "fips202.h"

/*************************************************
* Name:        hash_h
*
* Description: Hash function H used in KEM.
* Computes SHA3-256.
*
* Arguments:   - uint8_t *out:      pointer to output buffer (32 bytes)
* - const uint8_t *in: pointer to input buffer
* - size_t inlen:      length of input buffer in bytes
**************************************************/
static void hash_h(uint8_t *out, const uint8_t *in, size_t inlen) {
    sha3_256(out, in, inlen);
}

/*************************************************
* Name:        hash_g
*
* Description: Hash function G used in KEM.
* Computes SHAKE256 with an output of 2 * LORE_SYMBYTES bytes.
*
* Arguments:   - uint8_t *out:      pointer to output buffer (2 * LORE_SYMBYTES bytes)
* - const uint8_t *in: pointer to input buffer
* - size_t inlen:      length of input buffer in bytes
**************************************************/
static void hash_g(uint8_t *out, const uint8_t *in, size_t inlen) {
    shake256(out, 2 * LORE_SYMBYTES, in, inlen);
}

/*************************************************
* Name:        crypto_kem_keypair_derand
*
* Description: Generates public and private key for the CCA-secure KEM.
* (deterministic version)
*
* Arguments:   - unsigned char *pk:    pointer to output public key
* (an already allocated array of LORE_PUBLICKEYBYTES bytes)
* - unsigned char *sk:    pointer to output secret key
* (an already allocated array of LORE_SECRETKEYBYTES bytes)
* - const unsigned char *coins: pointer to input randomness
* (an already allocated array of 2 * LORE_SYMBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair_derand(unsigned char *pk, unsigned char *sk, const unsigned char *coins)
{
    // Call the IND-CPA keypair generation function
    indcpa_keypair_derand(pk, sk, coins);

    // Append public key, hash of public key, and z to the secret key
    unsigned char *sk_rest = sk + LORE_SECRETKEYBYTES;
    memcpy(sk_rest, pk, LORE_PUBLICKEYBYTES);
    sk_rest += LORE_PUBLICKEYBYTES;
    hash_h(sk_rest, pk, LORE_PUBLICKEYBYTES);
    sk_rest += LORE_SYMBYTES;
    memcpy(sk_rest, coins + LORE_SYMBYTES, LORE_SYMBYTES);

    return 0;
}

/*************************************************
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key for the CCA-secure KEM.
* (randomized version)
*
* Arguments:   - unsigned char *pk:    pointer to output public key
* (an already allocated array of LORE_PUBLICKEYBYTES bytes)
* - unsigned char *sk:    pointer to output secret key
* (an already allocated array of LORE_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
    unsigned char coins[2 * LORE_SYMBYTES];
    randombytes(coins, 2 * LORE_SYMBYTES);
    crypto_kem_keypair_derand(pk, sk, coins);
    return 0;
}

/*************************************************
* Name:        crypto_kem_enc_derand
*
* Description: Generates ciphertext and shared secret for a given public key.
* (deterministic version)
*
* Arguments:   - unsigned char *ct:    pointer to output ciphertext
* (an already allocated array of LORE_CIPHERTEXTBYTES bytes)
* - unsigned char *ss:    pointer to output shared secret
* (an already allocated array of LORE_SYMBYTES bytes)
* - const unsigned char *pk:    pointer to input public key
* (an already allocated array of LORE_PUBLICKEYBYTES bytes)
* - const unsigned char *coins: pointer to input randomness
* (an already allocated array of LORE_SYMBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_enc_derand(unsigned char *ct, unsigned char *ss, const unsigned char *pk, const unsigned char *coins)
{
    unsigned char kr[2 * LORE_SYMBYTES];
    unsigned char buf[2 * LORE_SYMBYTES];
    unsigned char mu[LORE_SYMBYTES];

    memcpy(mu, coins, LORE_SYMBYTES);

    hash_h(buf + LORE_SYMBYTES, pk, LORE_PUBLICKEYBYTES);
    memcpy(buf, mu, LORE_SYMBYTES);
    hash_g(kr, buf, 2 * LORE_SYMBYTES);

    indcpa_enc(ct, mu, pk, kr + LORE_SYMBYTES);

    memcpy(ss, kr, LORE_SYMBYTES);

    return 0;
}

/*************************************************
* Name:        crypto_kem_enc
*
* Description: Generates ciphertext and shared secret for a given public key.
* (randomized version)
*
* Arguments:   - unsigned char *ct:    pointer to output ciphertext
* (an already allocated array of LORE_CIPHERTEXTBYTES bytes)
* - unsigned char *ss:    pointer to output shared secret
* (an already allocated array of LORE_SYMBYTES bytes)
* - const unsigned char *pk:    pointer to input public key
* (an already allocated array of LORE_PUBLICKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
    unsigned char coins[LORE_SYMBYTES];
    randombytes(coins, LORE_SYMBYTES);
    crypto_kem_enc_derand(ct, ss, pk, coins);
    return 0;
}

/*************************************************
* Name:        crypto_kem_dec
*
* Description: Generates shared secret for a given ciphertext and private key.
*
* Arguments:   - unsigned char *ss:    pointer to output shared secret
* (an already allocated array of LORE_SYMBYTES bytes)
* - const unsigned char *ct:    pointer to input ciphertext
* (an already allocated array of LORE_CIPHERTEXTBYTES bytes)
* - const unsigned char *sk:    pointer to input secret key
* (an already allocated array of LORE_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{
    int fail;
    unsigned char mu[LORE_SYMBYTES];
    unsigned char kr[2 * LORE_SYMBYTES];
    unsigned char buf[2 * LORE_SYMBYTES];
    unsigned char ct_cmp[LORE_CIPHERTEXTBYTES];
    unsigned char ss_invalid[LORE_SYMBYTES];

    const unsigned char *sk_indcpa = sk;
    const unsigned char *pk = sk + LORE_SECRETKEYBYTES;
    const unsigned char *pkh = pk + LORE_PUBLICKEYBYTES;
    const unsigned char *z = pkh + LORE_SYMBYTES;

    indcpa_dec(mu, ct, sk_indcpa);

    memcpy(buf, mu, LORE_SYMBYTES);
    memcpy(buf + LORE_SYMBYTES, pkh, LORE_SYMBYTES);
    hash_g(kr, buf, 2 * LORE_SYMBYTES);

    indcpa_enc(ct_cmp, mu, pk, kr + LORE_SYMBYTES);

    fail = verify(ct, ct_cmp, LORE_CIPHERTEXTBYTES);

    rkprf(ss_invalid, z, ct);

    cmov(ss, kr, LORE_SYMBYTES, (unsigned char)(1 - fail));
    cmov(ss, ss_invalid, LORE_SYMBYTES, (unsigned char)fail);

    return 0;
}