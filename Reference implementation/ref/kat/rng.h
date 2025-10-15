#ifndef RNG_H
#define RNG_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
void
fprintBstr(FILE *fp, const char *S, const unsigned char *A, unsigned long long L);

int
FindMarker(FILE *fp, const char *marker);

int
ReadHex(FILE *fp, unsigned char *A, int Length, const char *str);

void
randombytes_init(unsigned char *entropy_input,
                 unsigned char *personalization_string,
                 int security_level);


void randombytes(uint8_t *out, size_t outlen);
uint32_t randombytes_uniform(uint32_t upper_bound);
#endif 