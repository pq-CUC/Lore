#include "symmetric.h"
#include "fips202.h"
#include <stdlib.h>

/*************************************************
* Name:        xof_absorb
*
* Description: Absorb step of the SHAKE-128 XOF.
* The domain-separation bytes are concatenated to the seed.
*
* Arguments:   - xof_state *s: pointer to (uninitialized) state
* - const uint8_t seed[LORE_SYMBYTES]: pointer to LORE_SYMBYTES input seed
* - uint8_t x:       additional byte of input
* - uint8_t y:       additional byte of input
**************************************************/
void xof_absorb(xof_state *s, const uint8_t seed[LORE_SYMBYTES], uint8_t x, uint8_t y) {
    uint8_t extseed[LORE_SYMBYTES + 2];

    for (size_t i = 0; i < LORE_SYMBYTES; ++i) {
        extseed[i] = seed[i];
    }
    extseed[LORE_SYMBYTES] = x;
    extseed[LORE_SYMBYTES + 1] = y;

    shake128_absorb_once(s, extseed, sizeof(extseed));
}

/*************************************************
* Name:        xof_squeezeblocks
*
* Description: Squeeze step of the SHAKE-128 XOF. Squeezes full blocks of
* 168 bytes each. Can be called multiple times to keep
* squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *out: pointer to output buffer
* - size_t outblocks: number of blocks to be squeezed
* - xof_state *s: pointer to input/output state
**************************************************/
void xof_squeezeblocks(uint8_t *out, size_t outblocks, xof_state *s) {
    shake128_squeezeblocks(out, outblocks, s);
}

/*************************************************
* Name:        prf
*
* Description: Usage of SHAKE-256 as a PRF, concatenating key and nonce.
*
* Arguments:   - uint8_t *out:       pointer to output buffer
* - size_t outlen:      number of requested output bytes
* - const uint8_t key[LORE_SYMBYTES]:  pointer to the key
* - uint8_t nonce:      single-byte nonce
**************************************************/
void prf(uint8_t *out, size_t outlen, const uint8_t key[LORE_SYMBYTES], uint8_t nonce) {
    uint8_t extkey[LORE_SYMBYTES + 1];

    for (size_t i = 0; i < LORE_SYMBYTES; ++i) {
        extkey[i] = key[i];
    }
    extkey[LORE_SYMBYTES] = nonce;

    shake256(out, outlen, extkey, sizeof(extkey));
}


/**
 * @brief Usage of SHAKE256 as a PRF for rejection sampling key generation.
 * Concatenates a secret key and a public input (ciphertext) to generate a pseudorandom key.
 *
 * @param[out] out      Pointer to the output buffer (LORE_SYMBYTES bytes).
 * @param[in]  key      Pointer to the secret key (LORE_SYMBYTES bytes).
 * @param[in]  input    Pointer to the public input, typically the ciphertext (LORE_CIPHERTEXTBYTES bytes).
 */
void lore_shake256_rkprf(unsigned char out[LORE_SYMBYTES], const unsigned char key[LORE_SYMBYTES], const unsigned char input[LORE_CIPHERTEXTBYTES])
{
  keccak_state s;

  shake256_init(&s);
  shake256_absorb(&s, key, LORE_SYMBYTES);
  shake256_absorb(&s, input, LORE_CIPHERTEXTBYTES);
  shake256_finalize(&s);
  shake256_squeeze(out, LORE_SYMBYTES, &s);
}