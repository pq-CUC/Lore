#include <string.h>
#include <ctype.h>
#include "rng.h"
#include "fips202.h"

static uint8_t      Key[32];
static uint8_t      V[16];
static int          reseed_counter = 0;

void
fprintBstr(FILE *fp, const char *S, const unsigned char *A, unsigned long long L)
{
    unsigned long long  i;
    fprintf(fp, "%s", S);
    for ( i=0; i<L; i++ )
        fprintf(fp, "%02X", A[i]);
    if ( L == 0 )
        fprintf(fp, "00");
    fprintf(fp, "\n");
}

int
FindMarker(FILE *fp, const char *marker)
{
    char    line[256];
    size_t  len;

    len = strlen(marker);
    if ( len > 254 )
        len = 254;

    while ( fgets(line, sizeof(line), fp) != NULL ) {
        if ( line[0] == '#' )
            continue;
        if ( strstr(line, marker) )
            return 1;
    }
    return 0;
}

int
ReadHex(FILE *fp, unsigned char *A, int Length, const char *str)
{
    int i, ch;
    unsigned char h;

    if ( Length == 0 ) {
        A[0] = 0x00;
        return 1;
    }

    if ( FindMarker(fp, str) == 0 )
        return 0;

    memset(A, 0x00, (size_t)Length);
    for ( i=0; i<Length*2; i++ ) {
        ch = fgetc(fp);
        while ( !isxdigit(ch) ) {
            if (ch == EOF) return 0;
            ch = fgetc(fp);
        }

        if (ch >= '0' && ch <= '9')
            h = (unsigned char)(ch - '0');
        else if (ch >= 'A' && ch <= 'F')
            h = (unsigned char)(ch - 'A' + 10);
        else if (ch >= 'a' && ch <= 'f')
            h = (unsigned char)(ch - 'a' + 10);
        else
            return 0;

        A[i/2] |= (unsigned char)(h << (4*(1-(i%2))));
    }

    return 1;
}


void
randombytes_init(unsigned char *entropy_input,
                 unsigned char *personalization_string,
                 int security_level)
{
    unsigned char   seed_material[48];

    (void)security_level;
    memcpy(seed_material, entropy_input, 48);
    if (personalization_string)
        for (int i=0; i<48; i++)
            seed_material[i] ^= personalization_string[i];
    memset(Key, 0x00, 32);
    memset(V, 0x00, 16);
    shake256(Key, 32, seed_material, 48);
    shake256(V, 16, seed_material, 48);
    reseed_counter = 1;
}

int
randombytes(unsigned char *x, unsigned long long xlen)
{
    unsigned char   block[16];
    int             i = 0;

    while ( xlen > 0 ) {
        for (int j=15; j>=0; j--) {
            if ( V[j] == 0xff )
                V[j] = 0x00;
            else {
                V[j]++;
                break;
            }
        }
        shake256(block, 16, V, 16);

        if ( xlen > 15 ) {
            memcpy(x+i, block, 16);
            i += 16;
            xlen -= 16;
        }
        else {
            memcpy(x+i, block, xlen);
            xlen = 0;
        }
    }
    reseed_counter++;

    return 1;
}

// ** Key fix: Add a deterministic implementation of randombytes_uniform **
uint32_t randombytes_uniform(uint32_t upper_bound)
{
  uint32_t r, min;

  if (upper_bound < 2) {
    return 0;
  }

  min = (uint32_t) -upper_bound % upper_bound;

  for (;;) {
    randombytes((unsigned char *)&r, sizeof(r));
    if (r >= min) {
      return r % upper_bound;
    }
  }
}