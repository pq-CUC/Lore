#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define ntt LORE_NAMESPACE(ntt)
void ntt(int16_t a[LORE_N]);

#define invntt_tomont LORE_NAMESPACE(invntt_tomont)
void invntt_tomont(int16_t a[LORE_N]);

#define fqmul LORE_NAMESPACE(fqmul)
int16_t fqmul(int16_t a, int16_t b);

#define basemul LORE_NAMESPACE(basemul)
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);

extern int16_t zetas[128];

#endif