#include "reduce.h"
#include "params.h"
#include <stdint.h>

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
* a*R^{-1} mod q, where R=2^16.
*
* Arguments:   - int64_t a: input integer to be reduced
*
* Returns:     montgomery_reduce(a) mod q.
**************************************************/
int16_t montgomery_reduce(int64_t a) {
    int32_t t;
    int16_t v;

    // QINV = -255, Q = 257
    v = (int16_t)(a * (-255));
    t = (int64_t)v * LORE_Q;
    t = (int32_t)(a - t);
    t >>= 16;
    return (int16_t)t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes a mod q
* in a way that is faster than a naive % operator.
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     a mod q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    // Use a 32-bit integer to prevent overflow in intermediate calculations
    int16_t t = a;
    
    // First reduction: t = (t & 0xFF) - (t >> 8)
    // This corresponds to X = r - k (since 2^8 ≡ -1 mod 257)
    t = (t & 0xFF) - (t >> 8);

    // Second reduction is commented out as per your request
     t = (t & 0xFF) - (t >> 8);
    
    return (int16_t)t;
}


/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial.
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
// Optimized version (unrolled 4 times)
void poly_reduce(poly *r) {
    for(int i=0; i<LORE_N; i+=4) {
        r->coeffs[i]   = barrett_reduce(r->coeffs[i]);
        r->coeffs[i+1] = barrett_reduce(r->coeffs[i+1]);
        r->coeffs[i+2] = barrett_reduce(r->coeffs[i+2]);
        r->coeffs[i+3] = barrett_reduce(r->coeffs[i+3]);
    }
}


