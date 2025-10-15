#include <stdio.h>
#include <string.h>
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "sampler.h"
#include "reduce.h"

/*************************************************
* Name:        pack_pk_ntt
*
* Description: Packs the public key into a byte array.
*
* Arguments:   - unsigned char *r:          pointer to the output byte array
* - const unsigned char *seedA:  pointer to the seed for matrix A
* - const poly_crt_vec *b_ntt:   pointer to the polynomial vector b in NTT domain
**************************************************/
static void pack_pk_ntt(unsigned char *r, const unsigned char *seedA, const poly_crt_vec *b_ntt) {
    // 1. Copy seedA
    memcpy(r, seedA, LORE_SYMBYTES);
    unsigned char *main_ptr = r + LORE_SYMBYTES;
    
    // 2. Prepare a temporary buffer large enough to collect all overflow bits
    unsigned char overflow_buf_all[LORE_K * LORE_N / 8 + 2] = {0};
    size_t total_overflow_bits = 0;

    // 3. Process K polynomials sequentially
    for (int i = 0; i < LORE_K; i++) {
        unsigned char temp_overflow_buf[LORE_N / 8 + 2] = {0};
        size_t bits_written = lore_poly_pack_q_split(main_ptr, temp_overflow_buf, &b_ntt->vec[i].q_poly);
        
        // Append the generated overflow bits to the total overflow buffer overflow_buf_all
        for (size_t k = 0; k < bits_written; k++) {
            if ((temp_overflow_buf[k / 8] >> (k % 8)) & 1) {
                size_t final_bit_pos = total_overflow_bits + k;
                overflow_buf_all[final_bit_pos / 8] |= (unsigned char)(1 << (final_bit_pos % 8));
            }
        }
        total_overflow_bits += bits_written;
        main_ptr += LORE_N; 
    }

    // 4. Write all collected overflow bits at once to the end of the main data
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    memcpy(main_ptr, overflow_buf_all, overflow_byte_len);

    // 5. Finally, write the T-data
#if LORE_R_BITS > 0
    unsigned char *t_ptr = main_ptr + overflow_byte_len;
    size_t t_poly_len = LORE_POLY_COMPRESSED_BYTES_T;
    for (int i = 0; i < LORE_K; i++) {
        uint8_t p_t_indices[LORE_N];
        for(int j=0; j<LORE_N; ++j) {
            p_t_indices[j] = (uint8_t)find_closest_r_idx(b_ntt->vec[i].t_poly.coeffs[j]);
        }
        pack_t_bits(t_ptr, p_t_indices);
        t_ptr += t_poly_len;
    }
#endif
}

/*************************************************
* Name:        unpack_pk_ntt
*
* Description: Unpacks the public key from a byte array.
*
* Arguments:   - poly_crt_vec *b_ntt:      pointer to the output polynomial vector b in NTT domain
* - unsigned char *seedA:     pointer to the output seed for matrix A
* - const unsigned char *pk:  pointer to the input byte array
**************************************************/
static void unpack_pk_ntt(poly_crt_vec *b_ntt, unsigned char *seedA, const unsigned char *pk) {
    // 1. Copy seedA
    memcpy(seedA, pk, LORE_SYMBYTES);
    const unsigned char *ptr = pk + LORE_SYMBYTES;
    
    // 2. [CORE] Scan the main data area to calculate the total number of overflow bits
    size_t main_data_len = LORE_K * LORE_N;
    size_t total_overflow_bits = 0;
    for(size_t i = 0; i < main_data_len; ++i) {
        if (ptr[i] == 0xFF) {
            total_overflow_bits++;
        }
    }
    
    // 3. Precisely locate the pointers for each data section
    const unsigned char *main_ptr = ptr;
    const unsigned char *overflow_ptr = ptr + main_data_len;

    int overflow_bit_pos = 0;

    // 4. Unpack Q-data
    for (int i = 0; i < LORE_K; i++) {
        lore_poly_unpack_q_split(&b_ntt->vec[i].q_poly, main_ptr + (i * LORE_N), overflow_ptr, &overflow_bit_pos);
    }

    // 5. Unpack T-data
#if LORE_R_BITS > 0
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    const unsigned char *t_ptr = overflow_ptr + overflow_byte_len; // <--- Correct starting address for T-data
    size_t t_poly_len = LORE_POLY_COMPRESSED_BYTES_T;
    for (int i = 0; i < LORE_K; i++) {
        uint8_t p_t_indices[LORE_N];
        unpack_t_bits(p_t_indices, t_ptr);
        t_ptr += t_poly_len;
        for (int j = 0; j < LORE_N; ++j) {
            b_ntt->vec[i].t_poly.coeffs[j] = (int16_t)(3 * p_t_indices[j] + 1);
        }
    }
#else
    for (int i = 0; i < LORE_K; i++) {
        for (int j = 0; j < LORE_N; j++) {
            b_ntt->vec[i].t_poly.coeffs[j] = 1;
        }
    }
#endif
}
/*************************************************
* Name:        pack_sk
*
* Description: Packs the secret key into a byte array.
*
* Arguments:   - unsigned char *r:             pointer to the output byte array
* - const poly_sparse *s_t_sparse: pointer to the sparse t-part of the secret key
* - const poly_crt_vec *s_crt_ntt: pointer to the CRT form of the secret key in NTT domain
**************************************************/
static void pack_sk(unsigned char *r, const poly_sparse *s_t_sparse, const poly_crt_vec *s_crt_ntt) {
    unsigned char *ptr = r;

    // 1. First, pack all n_coeffs
    for(int i = 0; i < LORE_K; i++) {
        memcpy(ptr, &s_t_sparse[i].n_coeffs, sizeof(int16_t));
        ptr += sizeof(int16_t);
    }

    // 2. Then, pack all valid pos and val pairs
    for (int i = 0; i < LORE_K; i++) {
        memcpy(ptr, s_t_sparse[i].pos, (size_t)s_t_sparse[i].n_coeffs * sizeof(uint8_t));
        ptr += (size_t)s_t_sparse[i].n_coeffs * sizeof(uint8_t);
        memcpy(ptr, s_t_sparse[i].val, (size_t)s_t_sparse[i].n_coeffs * sizeof(int8_t));
        ptr += (size_t)s_t_sparse[i].n_coeffs * sizeof(int8_t);
    }

    // 3. Finally, pack the q polynomials
    for (int i = 0; i < LORE_K; i++) {
        memcpy(ptr, s_crt_ntt->vec[i].q_poly.coeffs, LORE_POLY_BYTES);
        ptr += LORE_POLY_BYTES;
    }
}


/*************************************************
* Name:        unpack_sk
*
* Description: Unpacks the secret key from a byte array.
*
* Arguments:   - poly_crt_vec *s_crt_ntt_q:  pointer to the output CRT form of the secret key (q-part)
* - poly_sparse *s_t_sparse:    pointer to the output sparse t-part of the secret key
* - const unsigned char *sk:    pointer to the input byte array
**************************************************/
static void unpack_sk(poly_crt_vec *s_crt_ntt_q, poly_sparse *s_t_sparse, const unsigned char *sk) {
    const unsigned char *ptr = sk;

    // 1. First, unpack all n_coeffs
    for(int i = 0; i < LORE_K; i++) {
        memcpy(&s_t_sparse[i].n_coeffs, ptr, sizeof(int16_t));
        ptr += sizeof(int16_t);
    }

    // 2. Then, unpack all valid pos and val pairs
    for (int i = 0; i < LORE_K; i++) {
        memcpy(s_t_sparse[i].pos, ptr, (size_t)s_t_sparse[i].n_coeffs * sizeof(uint8_t));
        ptr += (size_t)s_t_sparse[i].n_coeffs * sizeof(uint8_t);
        memcpy(s_t_sparse[i].val, ptr, (size_t)s_t_sparse[i].n_coeffs * sizeof(int8_t));
        ptr += (size_t)s_t_sparse[i].n_coeffs * sizeof(int8_t);
    }

    // 3. Finally, unpack the q polynomials
    for (int i = 0; i < LORE_K; i++) {
        memcpy(s_crt_ntt_q->vec[i].q_poly.coeffs, ptr, LORE_POLY_BYTES);
        ptr += LORE_POLY_BYTES;
    }
}

/*************************************************
* Name:        pack_ciphertext_ntt
*
* Description: Packs the ciphertext into a byte array.
*
* Arguments:   - unsigned char *r:          pointer to the output byte array
* - const poly_crt_vec *cu_ntt:  pointer to the polynomial vector cu in NTT domain
* - const poly_crt *cv_crt:    pointer to the polynomial cv in CRT domain
**************************************************/
static void pack_ciphertext_ntt(unsigned char *r, const poly_crt_vec *cu_ntt, const poly_crt *cv_crt) {
    unsigned char *main_ptr = r;
    unsigned char overflow_buf_all[((LORE_K + 1) * LORE_N + 7) / 8 + 2] = {0};
    size_t total_overflow_bits = 0;
    
    // 1. Pack K cu polynomials
    for(int i = 0; i < LORE_K; ++i) {
        unsigned char temp_overflow_buf[LORE_N / 8 + 2] = {0};
        size_t bits_written = lore_poly_pack_q_split(main_ptr, temp_overflow_buf, &cu_ntt->vec[i].q_poly);
        
        for (size_t k = 0; k < bits_written; k++) {
            if ((temp_overflow_buf[k / 8] >> (k % 8)) & 1) {
                size_t final_bit_pos = total_overflow_bits + k;
                overflow_buf_all[final_bit_pos / 8] |= (unsigned char)(1 << (final_bit_pos % 8));
            }
        }
        total_overflow_bits += bits_written;
        main_ptr += LORE_N;
    }
    
    // 2. Pack 1 cv polynomial
    unsigned char temp_overflow_buf_cv[LORE_N / 8 + 2] = {0};
    size_t bits_written_cv = lore_poly_pack_q_split(main_ptr, temp_overflow_buf_cv, &cv_crt->q_poly);
    for (size_t k = 0; k < bits_written_cv; k++) {
        if ((temp_overflow_buf_cv[k / 8] >> (k % 8)) & 1) {
            size_t final_bit_pos = total_overflow_bits + k;
            overflow_buf_all[final_bit_pos / 8] |= (unsigned char)(1 << (final_bit_pos % 8));
        }
    }
    total_overflow_bits += bits_written_cv;
    main_ptr += LORE_N;

    // 3. Write the total overflow data
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    memcpy(main_ptr, overflow_buf_all, overflow_byte_len);

    // 4. Write T-data
#if LORE_R_BITS > 0
    unsigned char* t_ptr = main_ptr + overflow_byte_len;
    size_t t_poly_len = LORE_POLY_COMPRESSED_BYTES_T;
    for(int i = 0; i < LORE_K; ++i) {
        uint8_t t_indices[LORE_N];
        for(int j = 0; j < LORE_N; ++j) {
            t_indices[j] = (uint8_t)find_closest_r_idx(cu_ntt->vec[i].t_poly.coeffs[j]);
        }
        pack_t_bits(t_ptr, t_indices);
        t_ptr += t_poly_len;
    }
#endif
}

/*************************************************
* Name:        unpack_ciphertext_ntt
*
* Description: Unpacks the ciphertext from a byte array.
*
* Arguments:   - poly_crt_vec *cu_ntt:      pointer to the output polynomial vector cu in NTT domain
* - poly_crt *cv:             pointer to the output polynomial cv
* - const unsigned char *c:   pointer to the input byte array
**************************************************/
static void unpack_ciphertext_ntt(poly_crt_vec *cu_ntt, poly_crt *cv, const unsigned char *c) {
    // 1. The main data area of the ciphertext contains cu (K polynomials) and cv (1 polynomial)
    size_t main_data_len = (LORE_K + 1) * LORE_N;
    size_t total_overflow_bits = 0;
    for(size_t i = 0; i < main_data_len; ++i) {
        if (c[i] == 0xFF) {
            total_overflow_bits++;
        }
    }

    // 2. Determine the correct pointers for each data section
    const unsigned char *main_ptr = c;
    const unsigned char *overflow_ptr = c + main_data_len;

    int overflow_bit_pos = 0; // Used to track the bit position in the overflow stream

    // Unpack the q-part of K cu polynomials
    for (int i = 0; i < LORE_K; i++) {
        lore_poly_unpack_q_split(&cu_ntt->vec[i].q_poly, main_ptr + (i * LORE_N), overflow_ptr, &overflow_bit_pos);
    }
    // Unpack the q-part of the cv polynomial
    lore_poly_unpack_q_split(&cv->q_poly, main_ptr + (LORE_K * LORE_N), overflow_ptr, &overflow_bit_pos);

    // 3. Unpack t-poly from the correct position
#if LORE_R_BITS > 0
    size_t overflow_byte_len = (total_overflow_bits + 7) / 8;
    const unsigned char *t_ptr = overflow_ptr + overflow_byte_len; 
    size_t t_poly_len = LORE_POLY_COMPRESSED_BYTES_T; // This macro is defined in params.h
    for (int i = 0; i < LORE_K; i++) {
        uint8_t t_indices[LORE_N];
        unpack_t_bits(t_indices, t_ptr);
        t_ptr += t_poly_len;
        for(int j=0; j<LORE_N; ++j){
            cu_ntt->vec[i].t_poly.coeffs[j] = (int16_t)(3 * t_indices[j] + 1);
        }
    }
#else
    for(int i=0; i<LORE_K; ++i) {
        for(int j=0; j<LORE_N; ++j){
            cu_ntt->vec[i].t_poly.coeffs[j] = 1;
        }
    }
#endif

    // The t-part of cv is fixed, this logic remains unchanged
    const int16_t r0 = (LORE_T - 1) / 2;
    for (int i = 0; i < LORE_N; ++i) {
        cv->t_poly.coeffs[i] = r0;
    }
    poly_reduce(&cv->q_poly);
}


/*************************************************
* Name:        indcpa_keypair_derand
*
* Description: Deterministically generates IND-CPA public and private key.
* All randomness is derived from the input coins.
**************************************************/
void indcpa_keypair_derand(unsigned char *pk, unsigned char *sk, const unsigned char *coins)
{
    poly_crt_vec a_ntt[LORE_K];
    poly_crt_vec s_crt;
    poly_crt_vec b_final;
    poly_sparse s_t_sparse[LORE_K];

    unsigned char buf[2 * LORE_SYMBYTES];
    const unsigned char *seedA = buf;
    const unsigned char *noiseseed = buf + LORE_SYMBYTES;
    
    shake256(buf, 2 * LORE_SYMBYTES, coins, LORE_SYMBYTES);

    gen_matrix_ntt(a_ntt, seedA, 0);

    poly_getnoise(&s_crt, s_t_sparse, (unsigned char*)noiseseed, 0);

    poly_crt_vec_ntt(&s_crt);

    for (int i = 0; i < LORE_K; i++) {
    #if LORE_LEVEL == 1
        poly_crt_vec_pointwise_acc_montgomery_sparse(&b_final.vec[i], &a_ntt[i], &s_crt, s_t_sparse);
    #else
        poly_crt_vec_pointwise_acc_montgomery(&b_final.vec[i], &a_ntt[i], &s_crt);
    #endif
    }

    // (后续的舍入和打包逻辑保持不变)
    for(int i=0; i<LORE_K; ++i){
        poly delta_q_std;
        poly t_poly_orig = b_final.vec[i].t_poly;
        for(int j=0; j<LORE_N; ++j){
            int16_t x_t = t_poly_orig.coeffs[j];
            int16_t r_idx = find_closest_r_idx(x_t);
            int16_t r_val =(int16_t)( 3 * r_idx + 1);
            int16_t delta_val = r_val - x_t;
            while (delta_val > LORE_T / 2) delta_val -= LORE_T;
            while (delta_val <= -LORE_T / 2) delta_val += LORE_T;
            delta_q_std.coeffs[j] = delta_val;
        }

        poly_ntt(&delta_q_std);
        poly_add(&b_final.vec[i].q_poly, &b_final.vec[i].q_poly, &delta_q_std);

        for(int j=0; j<LORE_N; ++j){
            int16_t r_idx = find_closest_r_idx(t_poly_orig.coeffs[j]);
            b_final.vec[i].t_poly.coeffs[j] =(int16_t)( 3 * r_idx + 1);
        }
    }

    pack_pk_ntt(pk, seedA, &b_final);
    pack_sk(sk, s_t_sparse, &s_crt);
}

/*************************************************
* Name:        indcpa_keypair
*
* Description: Randomly generates IND-CPA public and private key.
* (This is now a wrapper around the deterministic version)
**************************************************/
void indcpa_keypair(unsigned char *pk, unsigned char *sk)
{
    unsigned char coins[LORE_SYMBYTES];
    randombytes(coins, LORE_SYMBYTES);
    indcpa_keypair_derand(pk, sk, coins);
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the IND-CPA secure
* public-key encryption scheme.
*
* Arguments:   - unsigned char *c:         pointer to output ciphertext
* - const unsigned char *m:   pointer to input message
* - const unsigned char *pk:  pointer to public key
* - const unsigned char *coins: pointer to randomness
**************************************************/
void indcpa_enc(unsigned char *c,
                const unsigned char *m,
                const unsigned char *pk,
                const unsigned char *coins) {
    poly_crt_vec at_ntt[LORE_K];
    poly_crt_vec sp_crt;
    poly_crt_vec b_ntt;
    poly_crt_vec cu_final;
    poly_crt cv_crt;
    poly e_pp; 
    poly_sparse sp_t_sparse[LORE_K];

    unsigned char seedA[LORE_SYMBYTES];

    unpack_pk_ntt(&b_ntt, seedA, pk);
    gen_matrix_ntt(at_ntt, seedA, 1); // Directly and efficiently generate the transpose matrix A^T

    poly_getnoise(&sp_crt, sp_t_sparse, (unsigned char*)coins, 0);

    poly_crt_vec_ntt(&sp_crt);

    for (int i = 0; i < LORE_K; i++) {
    #if LORE_LEVEL == 1
        // Level 1: Sparse multiplication is faster, keep as is
        poly_crt_vec_pointwise_acc_montgomery_sparse(&cu_final.vec[i], &at_ntt[i], &sp_crt, sp_t_sparse);
    #else
        // Level 2 & 3: The dense form of s' is available, call the dense multiplication function directly to avoid conversion
        poly_crt_vec_pointwise_acc_montgomery(&cu_final.vec[i], &at_ntt[i], &sp_crt);
    #endif
    }

    for(int i=0; i<LORE_K; ++i){
        poly delta_q_std;
        poly t_poly_orig = cu_final.vec[i].t_poly;
        for(int j=0; j<LORE_N; ++j){
            int16_t x_t = t_poly_orig.coeffs[j];
            int16_t r_idx = find_closest_r_idx(x_t);
            int16_t r_val = (int16_t)(3 * r_idx + 1);
            int16_t delta_val = r_val - x_t;
            while (delta_val > LORE_T / 2) delta_val -= LORE_T;
            while (delta_val <= -LORE_T / 2) delta_val += LORE_T;
            delta_q_std.coeffs[j] = delta_val;
        }
        
        poly_ntt(&delta_q_std);
        poly_add(&cu_final.vec[i].q_poly, &cu_final.vec[i].q_poly, &delta_q_std);

        for(int j=0; j<LORE_N; ++j){
             int16_t r_idx = find_closest_r_idx(t_poly_orig.coeffs[j]);
             cu_final.vec[i].t_poly.coeffs[j] = (int16_t)(3 * r_idx + 1);
        }
    }

    poly_crt_vec_pointwise_acc_montgomery_sparse(&cv_crt, (const poly_crt_vec*)&b_ntt, &sp_crt, sp_t_sparse);

    poly_invntt_tomont(&cv_crt.q_poly);
    for (int i = 0; i < LORE_N; ++i) {
        cv_crt.q_poly.coeffs[i] = montgomery_reduce((int64_t)cv_crt.q_poly.coeffs[i]);
    }
    poly_getnoise_uniform(&e_pp, LORE_T, coins, (unsigned char)LORE_K);
    for (int i = 0; i < LORE_N; ++i) {
        cv_crt.q_poly.coeffs[i] = barrett_reduce(cv_crt.q_poly.coeffs[i] + e_pp.coeffs[i]);
        int16_t t_coeff = (int16_t)((cv_crt.t_poly.coeffs[i] + e_pp.coeffs[i]) % LORE_T);
        if (t_coeff > LORE_T / 2) t_coeff -= LORE_T;
        if (t_coeff < -(LORE_T / 2)) t_coeff += LORE_T;
        cv_crt.t_poly.coeffs[i] = t_coeff;
    }

    poly msg_poly;
    poly_encode_msg(&msg_poly, m);
    poly_add_scaled_msg(&cv_crt, &msg_poly);

    // ======================== Start of fix ========================
    // Create a new poly_crt struct to hold the corrected cv data
    poly_crt cv_final_to_pack;

    // This is the core logic of lore_compress_cv_q
    const int16_t r0 = (LORE_T - 1) / 2;
    for (int i = 0; i < LORE_N; ++i) {
        int16_t x_q = cv_crt.q_poly.coeffs[i];
        int16_t x_t = cv_crt.t_poly.coeffs[i];
        int16_t delta = r0 - x_t;
        while (delta > LORE_T / 2) delta -= LORE_T;
        while (delta <= -LORE_T / 2) delta += LORE_T;
        // Store the corrected q-poly coefficients in the new struct
        cv_final_to_pack.q_poly.coeffs[i] = (int16_t)(x_q + delta + LORE_Q) % LORE_Q;
    }
    // The t-poly part is not passed to the packing function, but we copy it for completeness
    cv_final_to_pack.t_poly = cv_crt.t_poly;

    // Pass the corrected cv_final_to_pack (not the original cv_crt) to the packing function
    pack_ciphertext_ntt(c, &cu_final, &cv_final_to_pack);
    // ======================== End of fix ========================
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the IND-CPA secure
* public-key encryption scheme.
*
* Arguments:   - unsigned char *m:        pointer to output message
* - const unsigned char *c:  pointer to input ciphertext
* - const unsigned char *sk: pointer to secret key
**************************************************/
void indcpa_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk) {
    poly_crt_vec cu_ntt;
    poly_crt cv_crt, s_dot_cu;

    poly_crt_vec s_crt; // s_crt will be directly generated by the new unpack_sk
    poly_sparse s_t_sparse[LORE_K]; // New: for storing the unpacked sparse t-part


    unpack_ciphertext_ntt(&cu_ntt, &cv_crt, c);
    unpack_sk(&s_crt, s_t_sparse, sk);

    poly_crt_vec_pointwise_acc_montgomery_sparse(&s_dot_cu, (const poly_crt_vec*)&cu_ntt, &s_crt, s_t_sparse);

    poly_invntt_tomont(&s_dot_cu.q_poly);

    for (int i = 0; i < LORE_N; ++i) {
      s_dot_cu.q_poly.coeffs[i] = montgomery_reduce((int32_t)s_dot_cu.q_poly.coeffs[i]);
    }

    poly_crt R_crt;
    poly_crt_sub(&R_crt, &cv_crt, &s_dot_cu);

    poly_decode_msg_crt(m, &R_crt);
}