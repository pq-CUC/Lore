#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rng.h"
#include "api.h"

#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4

/*************************************************
* Name:        hex_str_to_bytes
*
* Description: Converts a hexadecimal string to a byte array.
*
* Arguments:   - unsigned char *dest: pointer to the destination byte array
* - const char *line:    pointer to the source hexadecimal string
* - int len:             the number of bytes to convert
*
* Returns 1 on success, 0 on failure.
**************************************************/
static int hex_str_to_bytes(unsigned char *dest, const char *line, int len) {
    const char *p = strstr(line, "= ");
    if (!p) return 0;
    p += 2; 

    for (int i = 0; i < len * 2; ++i) {
        char c1 = p[i];
        if (c1 == '\0' || c1 == '\n' || c1 == '\r') {
            return 0;
        }
        
        unsigned char byte_val;
        if (c1 >= '0' && c1 <= '9') byte_val = (unsigned char)(c1 - '0');
        else if (c1 >= 'a' && c1 <= 'f') byte_val = (unsigned char)(c1 - 'a' + 10);
        else if (c1 >= 'A' && c1 <= 'F') byte_val = (unsigned char)(c1 - 'A' + 10);
        else return 0;

        if (i % 2 == 0) {
            dest[i/2] = (byte_val << 4);
        } else {
            dest[i/2] |= byte_val;
        }
    }
    return 1;
}

/*************************************************
* Name:        main
*
* Description: Main function for the Known Answer Test (KAT) generation program.
* It creates request (.req) and response (.rsp) files for the PKE.
*
* Arguments:   - void
*
* Returns KAT_SUCCESS on success, or an error code on failure.
**************************************************/
int main() {
    char                fn_req[64], fn_rsp[64];
    FILE                *fp_req, *fp_rsp;
    unsigned char       seed[48];
    unsigned char       msg[CRYPTO_BYTES];
    unsigned char       entropy_input[48];
    int                 count;

    unsigned char       *ct, *pk, *sk;
    unsigned char       *pmsg;
    unsigned long long  mlen;

    sprintf(fn_req, "PQCpkeKAT_%s.req", CRYPTO_ALGNAME);
    if ((fp_req = fopen(fn_req, "w")) == NULL) {
        printf("Couldn't open <%s> for write\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }
    sprintf(fn_rsp, "PQCpkeKAT_%s.rsp", CRYPTO_ALGNAME);
    if ((fp_rsp = fopen(fn_rsp, "w")) == NULL) {
        printf("Couldn't open <%s> for write\n", fn_rsp);
        return KAT_FILE_OPEN_ERROR;
    }
    for (int i = 0; i < 48; i++) entropy_input[i] = (unsigned char)i;
    randombytes_init(entropy_input, NULL, 256);
    fprintf(fp_req, "# %s\n", CRYPTO_ALGNAME);
    for (int i = 0; i < 100; i++) {
        fprintf(fp_req, "count = %d\n", i);
        randombytes(seed, 48);
        fprintBstr(fp_req, "seed = ", seed, 48);
        fprintf(fp_req, "mlen = %d\n", CRYPTO_BYTES);
        randombytes(msg, CRYPTO_BYTES);
        fprintBstr(fp_req, "msg = ", msg, CRYPTO_BYTES);
        fprintf(fp_req, "\n");
    }
    fclose(fp_req);

    if ((fp_req = fopen(fn_req, "r")) == NULL) {
        printf("Couldn't open <%s> for read\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }
    fprintf(fp_rsp, "# %s\n", CRYPTO_ALGNAME);

    pk = (unsigned char *)calloc(CRYPTO_PUBLICKEYBYTES, sizeof(unsigned char));
    sk = (unsigned char *)calloc(CRYPTO_SECRETKEYBYTES, sizeof(unsigned char));
    ct = (unsigned char *)calloc(CRYPTO_CIPHERTEXTBYTES, sizeof(unsigned char));
    pmsg = (unsigned char *)calloc(CRYPTO_BYTES, sizeof(unsigned char));

    char line[512];
    while(fgets(line, sizeof(line), fp_req) && line[0] == '#');

    do {
        if (sscanf(line, "count = %d", &count) != 1) break;
        fprintf(fp_rsp, "count = %d\n", count);
        printf("DEBUG: Processing count %d...\n", count);

        if (!fgets(line, sizeof(line), fp_req) || !hex_str_to_bytes(seed, line, 48)) {
            fprintf(stderr, "ERROR: Failed to read seed for count %d\n", count); return KAT_DATA_ERROR;
        }
        fprintBstr(fp_rsp, "seed = ", seed, 48);
        randombytes_init(seed, NULL, 256);

        if (!fgets(line, sizeof(line), fp_req) || sscanf(line, "mlen = %llu", &mlen) != 1) {
            fprintf(stderr, "ERROR: Failed to read mlen for count %d\n", count); return KAT_DATA_ERROR;
        }

        if (!fgets(line, sizeof(line), fp_req) || !hex_str_to_bytes(msg, line, (int)mlen)) {
            fprintf(stderr, "ERROR: Failed to read msg for count %d\n", count); return KAT_DATA_ERROR;
        }
        
        if (crypto_pke_keypair(pk, sk) != 0) { printf("crypto_pke_keypair error\n"); return KAT_CRYPTO_FAILURE; }
        
        size_t pk_main_data_len = LORE_K * LORE_N;
        size_t pk_total_overflow_bits = 0;
        for(size_t i = 0; i < pk_main_data_len; ++i) {
            if (pk[LORE_SYMBYTES + i] == 0xFF) {
                pk_total_overflow_bits++;
            }
        }
        size_t pk_overflow_data_len = (pk_total_overflow_bits + 7) / 8;
        size_t pk_t_data_len = 0;
        #if LORE_R_BITS > 0
            pk_t_data_len = LORE_K * ((LORE_N * LORE_R_BITS + 7) / 8);
        #endif
        size_t actual_pk_len = LORE_SYMBYTES + pk_main_data_len + pk_overflow_data_len + pk_t_data_len;
        fprintBstr(fp_rsp, "pk = ", pk, actual_pk_len);

        fprintBstr(fp_rsp, "sk = ", sk, CRYPTO_SECRETKEYBYTES);
        
        if (crypto_pke_encrypt(ct, msg, mlen, pk) != 0) { printf("crypto_pke_encrypt error\n"); return KAT_CRYPTO_FAILURE; }
        
        size_t main_data_len = (LORE_K + 1) * LORE_N;
        size_t total_overflow_bits = 0;
        for(size_t i = 0; i < main_data_len; ++i) {
            if (ct[i] == 0xFF) {
                total_overflow_bits++;
            }
        }
        size_t overflow_data_len = (total_overflow_bits + 7) / 8;
        
        size_t t_data_len = 0;
        #if LORE_R_BITS > 0
            t_data_len = LORE_K * ((LORE_N * LORE_R_BITS + 7) / 8);
        #endif
        
        size_t actual_ct_len = main_data_len + overflow_data_len + t_data_len;
        
        fprintBstr(fp_rsp, "ct = ", ct, actual_ct_len);
        
        if (crypto_pke_decrypt(pmsg, &mlen, ct, sk) != 0) { printf("crypto_pke_decrypt error\n"); return KAT_CRYPTO_FAILURE; }
        
        if (memcmp(msg, pmsg, mlen)) { printf("Decrypted message does not match original message\n"); return KAT_CRYPTO_FAILURE; }
        fprintBstr(fp_rsp, "pmsg = ", pmsg, mlen);
        fprintf(fp_rsp, "\n");

    } while(fgets(line, sizeof(line), fp_req) && fgets(line, sizeof(line), fp_req));

    fclose(fp_req);
    fclose(fp_rsp);
    free(pk);
    free(sk);
    free(ct);
    free(pmsg);
    printf("KAT generation completed successfully.\n");
    return KAT_SUCCESS;
}