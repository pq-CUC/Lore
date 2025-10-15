#include <stdio.h>
#include <string.h>
#include "../api.h"
#include "../randombytes.h"

#define NTESTS 100

int main(void) {

    unsigned char pk[CRYPTO_PKE_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_PKE_SECRETKEYBYTES];
    unsigned char ct[CRYPTO_PKE_CIPHERTEXTBYTES];
    unsigned char msg[CRYPTO_BYTES];
    unsigned char pmsg[CRYPTO_BYTES];
    unsigned long long pmlen;
    int i;
    int fails = 0;

    printf("Starting tests for %s\n", CRYPTO_ALGNAME);

    for (i = 0; i < NTESTS; i++) {
        randombytes(msg, CRYPTO_BYTES);

        if (crypto_pke_keypair(pk, sk) != 0) {
            printf("ERROR in keypair on iteration %d\n", i);
            fails++;
            continue;
        }

        if (crypto_pke_encrypt(ct, msg, CRYPTO_BYTES, pk) != 0) {
            printf("ERROR in encrypt on iteration %d\n", i);
            fails++;
            continue;
        }

        if (crypto_pke_decrypt(pmsg, &pmlen, ct, sk) != 0) {
            printf("ERROR in decrypt on iteration %d\n", i);
            fails++;
            continue;
        }

        if (pmlen != CRYPTO_BYTES) {
            printf("ERROR: decrypted message length incorrect on iteration %d\n", i);
            fails++;
            continue;
        }
        if (memcmp(msg, pmsg, CRYPTO_BYTES)) {
            printf("ERROR: decrypted message does not match original on iteration %d\n", i);
            fails++;
            continue;
        }

        if((i+1) % 10 == 0) {
            printf("Finished test %d/%d\n", i + 1, NTESTS);
        }
    }

    if (fails) {
        printf("\nERROR: FAILED %d/%d tests!\n", fails, NTESTS);
    } else {
        printf("\nSUCCESS: PASSED all %d tests!\n", NTESTS);
    }

    return fails;
}