#include <stdint.h>
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include "fips202.h"
#include "poly.h"
#include "ntt.h"
#include "randombytes.h"
#include "reduce.h"
#include "sampler.h"
#include "symmetric.h"
#include "lore_luts.h" 
#include "toomcook.h"

/*************************************************
* Name:        poly_from_sparse
*
* Description: Converts a sparse polynomial to a dense polynomial.
*
* Arguments:   - poly *r:          pointer to the output dense polynomial
* - const poly_sparse *s: pointer to the input sparse polynomial
**************************************************/
void poly_from_sparse(poly *r, const poly_sparse *s) {
    memset(r->coeffs, 0, LORE_N * sizeof(int16_t));
    for (int i = 0; i < s->n_coeffs; ++i) {
        r->coeffs[s->pos[i]] = s->val[i];
    }
}


/*************************************************
* Name:        poly_add_modt
*
* Description: Adds two polynomials in the ring Z_t[X]/(X^N+1).
*
* Arguments:   - poly *r:       pointer to the output polynomial r = a + b
* - const poly *a: pointer to the first polynomial
* - const poly *b: pointer to the second polynomial
**************************************************/
void poly_add_modt(poly *r, const poly *a, const poly *b) {
    for(int j=0; j<LORE_N; ++j) {
        int16_t sum = a->coeffs[j] + b->coeffs[j];
        sum %= LORE_T;
        if (sum > LORE_T / 2) {
            sum -= LORE_T;
        }
        while (sum < -(LORE_T / 2)) {
            sum += LORE_T;
        }
        r->coeffs[j] = sum;
    }
}


/*************************************************
* Name:        poly_sparse_mul_modt
*
* Description: Multiplies a dense polynomial A with a sparse polynomial B,
* with the result being mod (t, x^n+1).
* (Based on a generalized implementation from SMAUG Fig. 1)
*
* Arguments:   - poly *r:             pointer to the output polynomial r = a_dense * b_sparse
* - const poly *a_dense:   pointer to the dense polynomial A
* - const poly_sparse *b_sparse: pointer to the sparse polynomial B
**************************************************/
void poly_sparse_mul_modt(poly *r, const poly *a_dense, const poly_sparse *b_sparse) {
    int32_t c_tmp[2 * LORE_N] = {0}; // Use a 32-bit temporary accumulator
    
    // 1. Accumulate products
    for (int i = 0; i < b_sparse->n_coeffs; i++) {
        uint8_t degree = b_sparse->pos[i];
        int8_t  value  = b_sparse->val[i];
        
        for (int j = 0; j < LORE_N; j++) {
            // Core operation: c[deg+j] += a[j] * value
            c_tmp[degree + j] += (int32_t)a_dense->coeffs[j] * value;
        }
    }
    
    // 2. Ring reduction (mod x^n+1) and mod t reduction
    for (int j = 0; j < LORE_N; j++) {
        // c[j] = c[j] - c[n+j]
        int32_t final_val = c_tmp[j] - c_tmp[LORE_N + j];

    #if LORE_LEVEL == 1
        // mod 3 (Level 1): use mod 3 sparse multiplication table
        r->coeffs[j] = LORE_MOD_T_SPARSE(final_val);
    #elif LORE_LEVEL == 2
        // mod 7 (Level 2): use mod 7 sparse multiplication table (activate new table)
        r->coeffs[j] = LORE_MOD_T_SPARSE(final_val);
    #else
        // LORE_LEVEL == 3: Manually inline Barrett reduction for extreme performance
        // Pre-calculated "magic number" M = floor(2^35 / 13)
        const int64_t M = 2643056797;
        
        // a * M >> 35 is an excellent approximation of a / 13
        int32_t q = (int32_t)(((int64_t)final_val * M) >> 35);
        
        // **Correction**: Rename local variable 'r' to 'reduced_val'
        int32_t reduced_val = final_val - q * 13;
        
        // Correction step 1: Ensure remainder is in the range [0, 12]
        if (reduced_val >= 13) reduced_val -= 13;
        if (reduced_val < 0) reduced_val += 13;

        // Correction step 2: Center the result to [-6, 6]
        if (reduced_val > 6) reduced_val -= 13;

        // **Correction**: Assign the correct result to the correct variable
        r->coeffs[j] = (int16_t)reduced_val;
    #endif
    }
}





/*************************************************
* Name:        inv
*
* Description: Computes the modular inverse using the Extended Euclidean Algorithm.
*
* Arguments:   - int16_t a: the integer
* - int16_t n: the modulus
*
* Returns:     the modular inverse of a modulo n, or -1 if not invertible
**************************************************/
static int16_t inv(int16_t a, int16_t n) {
    int16_t t = 0, newt = 1;
    int16_t r = n, newr = a;
    while (newr != 0) {
        int16_t quotient = r / newr;
        int16_t temp_t = t;
        t = newt;
        newt = (int16_t)(temp_t - quotient * newt);
        int16_t temp_r = r;
        r = newr;
        newr = (int16_t)(temp_r - quotient * newr);
    }
    if (r > 1) return -1; // not invertible
    if (t < 0) t = t + n;
    return t;
}

/*************************************************
* Name:        poly_crt_combine
*
* Description: Combine a polynomial in CRT form back to Z_tq.
* This version correctly handles centered coefficients.
*
* Arguments:   - poly *r:          pointer to the output polynomial
* - const poly_crt *pc: pointer to the input CRT polynomial
**************************************************/
void poly_crt_combine(poly *r, const poly_crt *pc) {
    const int16_t q_inv_t = inv(LORE_Q, LORE_T);

    for (int i = 0; i < LORE_N; ++i) {
        // a_q and a_t are centered coefficients
        int16_t a_q = pc->q_poly.coeffs[i];
        int16_t a_t = pc->t_poly.coeffs[i];
        
        // We need to find x such that:
        // x ≡ a_q (mod Q)
        // x ≡ a_t (mod T)

        // Step 1: Convert coefficients to the standard [0, n-1] range
        int32_t a_q_std = a_q;
        if (a_q_std < 0) a_q_std += LORE_Q;
        
        int32_t a_t_std = a_t;
        if (a_t_std < 0) a_t_std += LORE_T;
        
        // Step 2: Apply the standard Garner's CRT algorithm
        // x has the form x = a_q_std + k*Q
        // We need to solve for k
        // a_q_std + k*Q ≡ a_t_std (mod T)
        // k*Q ≡ a_t_std - a_q_std (mod T)
        // k ≡ (a_t_std - a_q_std) * q_inv_t (mod T)

        int32_t h = a_t_std - a_q_std;
        h %= LORE_T;
        if (h < 0) h += LORE_T;

        h = (h * q_inv_t) % LORE_T;
        
        // Step 3: Reconstruct the final result
        int32_t result = a_q_std + (int32_t)LORE_Q * h;
        r->coeffs[i] = (int16_t)result;
    }
}
/*************************************************
* Name:        poly_crt_decompose
*
* Description: Decompose a polynomial from Z_tq to centered CRT form.
*
* Arguments:   - poly_crt *pc:   pointer to the output CRT polynomial
* - const poly *p: pointer to the input polynomial
**************************************************/
void poly_crt_decompose(poly_crt *pc, const poly *p) {
    for (int i = 0; i < LORE_N; ++i) {
        pc->q_poly.coeffs[i] = barrett_reduce(p->coeffs[i]);

#if LORE_LEVEL == 1 || LORE_LEVEL == 2
        // Use lookup table defined for Level 1 or 2
        pc->t_poly.coeffs[i] = LORE_MOD_T_DECOMPOSE(p->coeffs[i]);
#else
        // Level 3 (mod 13) retains the original arithmetic operations
        int16_t t_coeff = p->coeffs[i] % LORE_T;
        if (t_coeff > LORE_T / 2) t_coeff -= LORE_T;
        if (t_coeff < -(LORE_T / 2)) t_coeff += LORE_T;
        pc->t_poly.coeffs[i] = t_coeff;
#endif
    }
}

/*************************************************
* Name:        poly_crt_add
*
* Description: Adds two CRT polynomials.
*
* Arguments:   - poly_crt *r:       pointer to the output CRT polynomial
* - const poly_crt *a: pointer to the first input CRT polynomial
* - const poly_crt *b: pointer to the second input CRT polynomial
**************************************************/
void poly_crt_add(poly_crt *r, const poly_crt *a, const poly_crt *b) {
    for (int i = 0; i < LORE_N; ++i) {
        r->q_poly.coeffs[i] = (int16_t)((a->q_poly.coeffs[i] + b->q_poly.coeffs[i]) % LORE_Q);
        r->t_poly.coeffs[i] = (int16_t)((a->t_poly.coeffs[i] + b->t_poly.coeffs[i]) % LORE_T);
    }
}

/*************************************************
* Name:        poly_crt_sub
*
* Description: Subtracts two CRT polynomials.
*
* Arguments:   - poly_crt *r:       pointer to the output CRT polynomial
* - const poly_crt *a: pointer to the first input CRT polynomial
* - const poly_crt *b: pointer to the second input CRT polynomial
**************************************************/
void poly_crt_sub(poly_crt *r, const poly_crt *a, const poly_crt *b) {
    for (int i = 0; i < LORE_N; ++i) {
        r->q_poly.coeffs[i] = barrett_reduce(a->q_poly.coeffs[i] - b->q_poly.coeffs[i]);

        int16_t t_coeff = (int16_t)((a->t_poly.coeffs[i] - b->t_poly.coeffs[i]) % LORE_T);
        if (t_coeff > LORE_T / 2) t_coeff -= LORE_T;
        if (t_coeff < -(LORE_T / 2)) t_coeff += LORE_T;
        r->t_poly.coeffs[i] = t_coeff;
    }
}

/*************************************************
* Name:        poly_add
*
* Description: Adds two polynomials.
*
* Arguments:   - poly *r:       pointer to the output polynomial
* - const poly *a: pointer to the first input polynomial
* - const poly *b: pointer to the second input polynomial
**************************************************/
 void poly_add(poly *r, const poly *a, const poly *b) {
     for(int i=0; i<LORE_N; ++i) {
         r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
     }
     for(int i=0; i<LORE_N; ++i) {
         r->coeffs[i] = barrett_reduce(r->coeffs[i]);
     }
 }


/*************************************************
* Name:        poly_frombytes
*
* Description: Deserializes a polynomial from a byte array.
*
* Arguments:   - poly *r:            pointer to the output polynomial
* - const unsigned char *a: pointer to the input byte array
**************************************************/
void poly_frombytes(poly *r, const unsigned char *a) {
    for (int i = 0; i < LORE_N / 2; i++) {
        r->coeffs[2 * i] = ((a[3 * i + 0] >> 0) | ((int16_t)a[3 * i + 1] << 8)) & 0xFFF;
        r->coeffs[2 * i + 1] = ((a[3 * i + 1] >> 4) | ((int16_t)a[3 * i + 2] << 4)) & 0xFFF;
    }
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serializes a polynomial into a byte array.
*
* Arguments:   - unsigned char *r: pointer to the output byte array
* - const poly *p:    pointer to the input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, const poly *p) {
    for (int i = 0; i < LORE_N / 2; i++) {
        r[3 * i + 0] = (unsigned char)(p->coeffs[2 * i + 0] >> 0);
        r[3 * i + 1] = (unsigned char)(((p->coeffs[2 * i + 0] >> 8) & 0xFF) | (p->coeffs[2 * i + 1] << 4));
        r[3 * i + 2] = (unsigned char)(p->coeffs[2 * i + 1] >> 4);
    }
}

/*************************************************
* Name:        poly_getnoise
*
* Description: Samples a noise polynomial with a fixed weight.
*
* Arguments:   - poly_crt_vec *r_crt_vec:   pointer to the output CRT polynomial vector
* - poly_sparse *r_sparse_vec: pointer to the output sparse polynomial vector
* - unsigned char *seed:       pointer to the seed
* - unsigned char nonce:       a domain-separation nonce
**************************************************/
void poly_getnoise(poly_crt_vec *r_crt_vec, poly_sparse *r_sparse_vec, unsigned char *seed, unsigned char nonce)
{
  sample_fixed_weight(r_crt_vec, r_sparse_vec, seed, nonce);
}

/*************************************************
* Name:        poly_getnoise_uniform
*
* Description: Sample a polynomial with coefficients uniformly random from
* U[-(t-1)/2, (t-1)/2].
*
* Arguments:   - poly *r:                pointer to output polynomial
* - uint16_t t:             the parameter t
* - const unsigned char *seed: pointer to input seed
* - unsigned char nonce:      a domain-separation nonce
**************************************************/
void poly_getnoise_uniform(poly *r, uint16_t t, const unsigned char *seed, unsigned char nonce)
{
    // --- Optimization Core: Streamlined Sampling ---
    // 1. Generate enough random bytes at once.
    // To handle rejection sampling, we generate slightly more bytes than needed. LORE_N * 2 = 512 bytes is more than enough.
    uint8_t buf[LORE_N * 2];
    prf(buf, sizeof(buf), seed, nonce);

    int ctr = 0; // counter for successfully sampled coefficients
    unsigned int buf_pos = 0; // current position in the buffer

    // 2. Efficient rejection sampling threshold.
    // To avoid modular bias, we only accept random bytes that fall within the largest multiple of t.
    // e.g., for t=7, the threshold limit = floor(256/7)*7 = 36*7 = 252.
    // Any byte less than 252 can be taken modulo 7 without bias.
    uint16_t limit = (uint16_t)((0x100 / t) * t);

    // 3. Loop in memory, extracting bytes from the buffer and performing sampling.
    while(ctr < LORE_N) {
        // If the buffer is exhausted (though nearly impossible in this case), abort as a precaution.
        if (buf_pos >= sizeof(buf)) {
            break; 
        }

        uint8_t val = buf[buf_pos++];
        
        // 4. Perform rejection sampling.
        // The acceptance rate is very high (e.g., 252/256 ≈ 98.4% for t=7), much faster than external calls.
        if (val < limit) {
            // Sampling successful, center and assign
        #if LORE_LEVEL == 1 || LORE_LEVEL == 2
            // Use lookup table defined for Level 1 or 2
            // (the table generated by our python script is already centered)
            r->coeffs[ctr++] = LORE_MOD_T_U8_NOISE(val);
        #else
            // Level 3 (mod 13) retains the original arithmetic operations
            int16_t bound = (int16_t)((t - 1) / 2);
            // Sampling successful, center and assign
            r->coeffs[ctr++] = (int16_t)((val % t) - bound);
        #endif
        }
    }
}

/*************************************************
* Name:        poly_ntt
*
* Description: Applies Number Theoretic Transform (NTT) to a polynomial.
*
* Arguments:   - poly *r: pointer to the polynomial
**************************************************/
void poly_ntt(poly *r) {
    ntt(r->coeffs);
    poly_reduce(r);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Applies inverse NTT to a polynomial.
*
* Arguments:   - poly *r: pointer to the polynomial
**************************************************/
void poly_invntt_tomont(poly *r) {
    invntt_tomont(r->coeffs);
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of two polynomials in the NTT domain.
*
* Arguments:   - poly *r:       pointer to the output polynomial
* - const poly *a: pointer to the first input polynomial
* - const poly *b: pointer to the second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *r, const poly *a, const poly *b) {
    unsigned int i;
    for(i=0; i<LORE_N/4; i++) {
      basemul(&r->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas[64+i]);
      basemul(&r->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2], -zetas[64+i]);
    }
}

/*************************************************
* Name:        pack_t_bits
*
* Description: Packs the t-part indices into a byte array.
*
* Arguments:   - unsigned char *r:      pointer to the output byte array
* - const uint8_t *t_indices: pointer to the input t_indices
**************************************************/
#if LORE_R_BITS > 0
void pack_t_bits(unsigned char *r, const uint8_t *t_indices) {
    memset(r, 0, LORE_POLY_COMPRESSED_BYTES_T);
    for (size_t i = 0; i < LORE_N; i++) {
        for (size_t j = 0; j < LORE_R_BITS; j++) {
            if (((t_indices[i] >> j) & 1)) {
                r[(i * LORE_R_BITS + j) / 8] |= (unsigned char)(1 << ((i * LORE_R_BITS + j) % 8));
            }
        }
    }
}

/*************************************************
* Name:        unpack_t_bits
*
* Description: Unpacks the t-part indices from a byte array.
*
* Arguments:   - uint8_t *t_indices:    pointer to the output t_indices
* - const unsigned char *r: pointer to the input byte array
**************************************************/
void unpack_t_bits(uint8_t *t_indices, const unsigned char *r) {
    memset(t_indices, 0, LORE_N);
    for (size_t i = 0; i < LORE_N; i++) {
        for (size_t j = 0; j < LORE_R_BITS; j++) {
            if (((r[(i * LORE_R_BITS + j) / 8] >> ((i * LORE_R_BITS + j) % 8)) & 1)) {
               t_indices[i] |= (uint8_t)(1 << j);
            }
        }
    }
}
#endif

/*************************************************
* Name:        find_closest_r_idx (OPTIMIZED VERSION)
*
* Description: Finds the index of the closest value in the set R = {3i+1}
* to a given coefficient x_t mod t.
* This version uses a tiny, cache-friendly lookup table.
*
* Arguments:   - int16_t x_t: the input coefficient
*
* Returns:     the index of the closest value in R
**************************************************/
int16_t find_closest_r_idx(int16_t x_t) {
    #if LORE_R_BITS == 0
        // Level 1: |R|=1, index is always 0.
        (void)x_t;
        return 0;

    #else
        // CRITICAL STEP: First, reduce x_t to the centered range [-T/2, T/2],
        // mimicking the implicit behavior of the original loop's distance calculation.
        while (x_t > LORE_T / 2) x_t -= LORE_T;
        while (x_t < -(LORE_T / 2)) x_t += LORE_T;

        #if LORE_LEVEL == 2
            // Level 2: T=7, R={1, 4}. Centered x_t is in [-3, 3].
            // Access: lut[x_t + 3]
            static const int16_t lut_mod7[LORE_T] = {
                // x_t:     -3, -2, -1,  0,  1,  2,  3
                // R_idx:
                           1,  1,  0,  0,  0,  0,  1
            };
            return lut_mod7[x_t + 3];

        #elif LORE_LEVEL == 3
            // Level 3: T=13, R={1, 4, 7, 10}. Centered x_t is in [-6, 6].
            // Access: lut[x_t + 6]
            static const int16_t lut_mod13[LORE_T] = {
                // x_t:     -6, -5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6
                // R_idx:
                           2,  2,  3,  3,  0,  0,  0,  0,  0,  1,  1,  2,  2
            };
            return lut_mod13[x_t + 6];
        #endif
    #endif
}




#if LORE_LEVEL == 3 // t=13
    static const int16_t Q_INV_T = 4; // inv(257, 13)
#elif LORE_LEVEL == 2 // t=7
    static const int16_t Q_INV_T = 3; // inv(257, 7)
#else // LORE_LEVEL == 1, t=3
    static const int16_t Q_INV_T = 2; // inv(257, 3)
#endif

/*************************************************
* Name:        poly_decode_msg_crt
*
* Description: Decodes a message from a CRT polynomial.
*
* Arguments:   - unsigned char *msg:   pointer to the output message
* - const poly_crt *r_crt: pointer to the input CRT polynomial
**************************************************/
void poly_decode_msg_crt(unsigned char *msg, const poly_crt *r_crt) {
    const int32_t TQ = LORE_T * LORE_Q;
    const int32_t TQ_HALF = TQ / 2;
    const int32_t TQ_QUARTER = TQ / 4;

    memset(msg, 0, LORE_SYMBYTES);

    for (int i = 0; i < LORE_N; ++i) {
        int16_t a_q = r_crt->q_poly.coeffs[i];
        int16_t a_t = r_crt->t_poly.coeffs[i];

        // Convert centered coefficients back to the standard range [0, m-1]
        int32_t a_q_std = (a_q < 0) ? (a_q + LORE_Q) : a_q;
        int32_t a_t_std = (a_t < 0) ? (a_t + LORE_T) : a_t;

        // This is the core of Garner's algorithm for reconstruction
        // h = ((a_t - a_q) * q_inv_t) mod t
        int32_t h = a_t_std - a_q_std;
        h = h * Q_INV_T;

        // **** This is the key fix ****
        // The '%' operator in C can yield a negative result for negative numbers,
        // so it must be corrected to ensure h is in the range [0, t-1]
        h = (h % LORE_T + LORE_T) % LORE_T;

        // Reconstruct the full coefficient: R = a_q + h * q
        int32_t val = a_q_std + (int32_t)LORE_Q * h;

        // Center val to [-TQ/2, TQ/2)
        if (val > TQ_HALF) {
            val -= TQ;
        }

        // Decode using the same logic as the original poly_decode_msg
        if (abs(val) > TQ_QUARTER) {
            msg[i / 8] |= (unsigned char)(1 << (i % 8));
        }
    }
}


/*************************************************
* Name:        poly_encode_msg (Algorithm 9)
*
* Description: Encodes a message into a polynomial.
*
* Arguments:   - poly *r:                pointer to the output polynomial
* - const unsigned char *msg: pointer to the input message
**************************************************/
void poly_encode_msg(poly *r, const unsigned char *msg) {
    // Calculate the encoding value based on the paper's formula (tq-1)/2
    const int32_t val = (LORE_T * LORE_Q - 1) / 2;

    for (int i = 0; i < LORE_N / 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            uint16_t mask = 1 << j;
            if (msg[i] & mask) {
                r->coeffs[8*i+j] = (int16_t)val;
            } else {
                r->coeffs[8*i+j] = 0;
            }
        }
    }
}

/*************************************************
* Name:        poly_decode_msg (Algorithm 11)
*
* Description: Decodes a message from a polynomial.
*
* Arguments:   - unsigned char *msg: pointer to the output message
* - const poly *r:      pointer to the input polynomial
**************************************************/
void poly_decode_msg(unsigned char *msg, const poly *r) {
    const int32_t TQ = LORE_T * LORE_Q;
    const int32_t threshold = TQ / 4;

    memset(msg, 0, LORE_SYMBYTES);

    for (int i = 0; i < LORE_N; ++i) {
        // Center the coefficient from [0, TQ-1] to [-TQ/2, TQ/2)
        int32_t val = r->coeffs[i];
        if (val > TQ / 2) {
            val -= TQ;
        }

        // If the coefficient is closer to 0 than to +/- TQ/2, the bit is 0.
        // If it is closer to +/- TQ/2, the bit is 1.
        if (abs(val) > threshold) {
            msg[i / 8] |= (unsigned char)(1 << (i % 8));
        }
    }
}
/*************************************************
* Name:        poly_add_scaled_msg
*
* Description: Adds a message to a CRT polynomial. (Optimized version)
*
* Arguments:   - poly_crt *r:          pointer to the CRT polynomial
* - const poly *msg_poly: pointer to the message polynomial
**************************************************/
void poly_add_scaled_msg(poly_crt *r, const poly *msg_poly) {
    for(int i = 0; i < LORE_N; ++i) {
        // --- q-part addition (unchanged) ---
        r->q_poly.coeffs[i] = barrett_reduce(r->q_poly.coeffs[i] + msg_poly->coeffs[i]);

        // --- t-part addition (now highly optimized) ---
        int16_t t_msg_coeff = msg_poly->coeffs[i] % LORE_T;
        int32_t t_sum = r->t_poly.coeffs[i] + t_msg_coeff; // Use int32_t for sum

#if LORE_LEVEL == 1 || LORE_LEVEL == 2
        // For t=3 and t=7, use the extremely fast lookup table.
        // The LORE_MOD_T_SPARSE macro handles both modulo and centering.
        r->t_poly.coeffs[i] = LORE_MOD_T_SPARSE(t_sum);

#else // LORE_LEVEL == 3 (t=13)
        // For t=13, use the highly optimized Barrett-like reduction
        // adapted from poly_sparse_mul_modt function.

        // Pre-calculated "magic number" M = floor(2^35 / 13)
        const int64_t M = 2643056797;

        // a * M >> 35 is an excellent approximation of a / 13
        int32_t q = (int32_t)(((int64_t)t_sum * M) >> 35);
        int16_t reduced_val = (int16_t)(t_sum - q * 13);

        // Correction step 1: Ensure remainder is in the range [0, 12]
        if (reduced_val >= 13) reduced_val -= 13;
        if (reduced_val < 0) reduced_val += 13;

        // Correction step 2: Center the result to [-6, 6]
        if (reduced_val > 6) reduced_val -= 13;

        r->t_poly.coeffs[i] = reduced_val;
#endif
    }
}

/*************************************************
* Name:        rej_uniform_q
*
* Description: Rejection sampling for q-part coefficients.
* Consumes bytes from a buffer and generates uniform integers mod Q.
*
* Arguments:   - int16_t *r:           pointer to output buffer
* - unsigned int len:     requested number of integers
* - const uint8_t *buf:   pointer to input buffer
* - unsigned int buflen:  length of input buffer in bytes
*
* Returns:     number of sampled integers (at most len)
**************************************************/
unsigned int rej_uniform_q(int16_t *r,
                           unsigned int len,
                           const uint8_t *buf,
                           unsigned int buflen)
{
  unsigned int ctr = 0, pos = 0;
  uint16_t val;

  while(ctr < len && pos + 2 <= buflen) {
    val = (uint16_t)(buf[pos] | ((uint16_t)buf[pos+1] << 8));
    pos += 2;

    if (val < 65535) {
      r[ctr++] = barrett_reduce((int16_t)val);
    }
  }
  return ctr;
}

/*************************************************
* Name:        rej_uniform_t
*
* Description: Rejection sampling for t-part coefficients.
* Consumes bytes from a buffer and generates uniform integers mod T.
*
* Arguments:   - int16_t *r:           pointer to output buffer
* - unsigned int len:     requested number of integers
* - const uint8_t *buf:   pointer to input buffer
* - unsigned int buflen:  length of input buffer in bytes
*
* Returns:     number of sampled integers (at most len)
**************************************************/
unsigned int rej_uniform_t(int16_t *r,
                           unsigned int len,
                           const uint8_t *buf,
                           unsigned int buflen)
{
  unsigned int ctr = 0, pos = 0;
  const uint16_t t_limit = (0x100 / LORE_T) * LORE_T;

  while(ctr < len && pos < buflen) {
    if (buf[pos] < t_limit) {
      int16_t t_coeff = buf[pos] % LORE_T;
      if (t_coeff > LORE_T / 2) {
          t_coeff -= LORE_T;
      }
      r[ctr++] = t_coeff;
    }
    pos++;
  }
  return ctr;
}
/*************************************************
* Name:        poly_mul_modt
*
* Description: Karatsuba multiplication of two polynomials in Z_t[X]/(X^N+1).
*
* Arguments:   - poly *r:       pointer to the output polynomial
* - const poly *a: pointer to the first input polynomial
* - const poly *b: pointer to the second input polynomial
**************************************************/
void poly_mul_modt(poly *r, const poly *a, const poly *b) {
    // Step 1: Use the most efficient Toom-Cook algorithm for the core multiplication
    // This is an accumulation operation (res += a*b), so we clear the result first
    int16_t res_full[LORE_N];
    memset(res_full, 0, LORE_N * sizeof(int16_t));
    poly_mul_acc(a->coeffs, b->coeffs, res_full);

    // Step 2: Use the most efficient, level-optimized method for mod t reduction
    for (int i = 0; i < LORE_N; i++) {
        int32_t final_val = res_full[i]; // Use a 32-bit integer to match the input range of the lookup table

#if LORE_LEVEL == 1
        // Level 1 (mod 3): use the extremely fast lookup table
        r->coeffs[i] = LORE_MOD_T_SPARSE(final_val);
#elif LORE_LEVEL == 2
        // Level 2 (mod 7): use the extremely fast lookup table
        r->coeffs[i] = LORE_MOD_T_SPARSE(final_val);
#else
        // Level 3 (mod 13): use efficient arithmetic operations for centered reduction
        // (because no lookup table was designed for Level 3)
        int16_t reduced_val = (int16_t)(final_val % LORE_T);
        if (reduced_val > LORE_T / 2) {
            reduced_val -= LORE_T;
        } else if (reduced_val < -(LORE_T / 2)) {
            reduced_val += LORE_T;
        }
        r->coeffs[i] = reduced_val;
#endif
    }
}

/*************************************************
* Name:        print_poly
*
* Description: Prints the coefficients of a polynomial for debugging.
*
* Arguments:   - const char* name: the name of the polynomial
* - const poly* p:      the polynomial to print
**************************************************/
void print_poly(const char* name, const poly* p) {
    printf("--- Poly: %s ---\n", name);
    for (int i = 0; i < LORE_N; ++i) {
        printf("%5d ", p->coeffs[i]);
        if ((i + 1) % 16 == 0) {
            printf("\n");
        }
    }
    printf("-----------------------------------\n");
}

/*************************************************
* Name:        print_poly_crt
*
* Description: Prints the coefficients of a CRT polynomial for debugging.
*
* Arguments:   - const char* name:   the name of the polynomial
* - const poly_crt* p_crt: the CRT polynomial to print
**************************************************/
void print_poly_crt(const char* name, const poly_crt* p_crt) {
    printf("--- Poly CRT: %s ---\n", name);
    printf("q_poly:\n");
    for (int i = 0; i < LORE_N; ++i) {
        printf("%5d ", p_crt->q_poly.coeffs[i]);
        if ((i + 1) % 16 == 0) {
            printf("\n");
        }
    }
    printf("\nt_poly:\n");
    for (int i = 0; i < LORE_N; ++i) {
        printf("%5d ", p_crt->t_poly.coeffs[i]);
        if ((i + 1) % 16 == 0) {
            printf("\n");
        }
    }
    printf("-----------------------------------\n");
}


/*************************************************
* Name:        lore_poly_pack_q_split
*
* Description: Splits the q-part of a polynomial into a main buffer and an overflow buffer.
*
* Arguments:   - unsigned char main_buf[LORE_N]:      output, the fixed-size main data buffer
* - unsigned char *overflow_buf:       output, buffer for overflow bits
* - const poly *p:                     input, the polynomial to pack
*
* Returns:     the number of overflow bits generated
**************************************************/
size_t lore_poly_pack_q_split(unsigned char main_buf[LORE_N], unsigned char *overflow_buf, const poly *p) {
    int total_overflow_bits = 0;

    for (int j = 0; j < LORE_N; ++j) {
        int16_t coeff = p->coeffs[j];
        if (coeff < 0 || coeff >= LORE_Q) {
           coeff = barrett_reduce(coeff);
           if(coeff < 0) coeff += LORE_Q;
        }

        if (coeff < 255) {
            main_buf[j] = (unsigned char)coeff;
        } else {
            main_buf[j] = 0xFF; // Mark as overflow

            // Ensure the target byte is clean before writing the first overflow bit
            if ((total_overflow_bits % 8) == 0) {
                overflow_buf[total_overflow_bits / 8] = 0;
            }
            
            if (coeff == 256) {
                // If the coefficient is 256, set the corresponding bit to 1
                overflow_buf[total_overflow_bits / 8] |= (unsigned char)(1 << (total_overflow_bits % 8));
            }
            // If the coefficient is 255, we do nothing; the bit remains 0
            total_overflow_bits++;
        }
    }
    return (size_t)total_overflow_bits;
}


/*************************************************
* Name:        lore_poly_unpack_q_split
*
* Description: Reconstructs a polynomial from a main buffer and an overflow buffer.
*
* Arguments:   - poly *p:                         output, the reconstructed polynomial
* - const unsigned char main_buf[LORE_N]: input, the main data buffer
* - const unsigned char *overflow_buf:   input, the stream of overflow bits
* - int *overflow_bit_pos:           input/output, the current bit position in the overflow stream
**************************************************/
void lore_poly_unpack_q_split(poly *p, const unsigned char main_buf[LORE_N], const unsigned char *overflow_buf, int *overflow_bit_pos) {
    for (int j = 0; j < LORE_N; ++j) {
        unsigned char byte = main_buf[j];
        if (byte == 0xFF) {
            int bit = (overflow_buf[*overflow_bit_pos / 8] >> (*overflow_bit_pos % 8)) & 1;
            (*overflow_bit_pos)++;
            p->coeffs[j] = (int16_t)(255 + bit);
        } else {
            p->coeffs[j] = byte;
        }
    }
}