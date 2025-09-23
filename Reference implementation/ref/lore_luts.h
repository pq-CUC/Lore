#ifndef LORE_LUTS_H
#define LORE_LUTS_H

// params.h must be included before checking LORE_LEVEL.
#include "params.h" 
#include <stdint.h>

#if LORE_LEVEL == 1
    #include "lut_mod3.h"

    #define LORE_MOD_T_SPARSE(x) (LORE_LUT_SPARSE_MOD3[(x) + 160])
    
    #define LORE_MOD_T_DECOMPOSE(x) (LORE_LUT_DECOMPOSE_MOD3[(x) + 1])
    
    #define LORE_MOD_T_U8_NOISE(x) (LORE_LUT_U8_NOISE_MOD3[x])

#elif LORE_LEVEL == 2
    #include "lut_mod7.h"

    // Offset +1224
    #define LORE_MOD_T_SPARSE(x) (LORE_LUT_SPARSE_MOD7[(x) + 1224])
    
    // Offset +2
    #define LORE_MOD_T_DECOMPOSE(x) (LORE_LUT_DECOMPOSE_MOD7[(x) + 2])
    
    #define LORE_MOD_T_U8_NOISE(x) (LORE_LUT_U8_NOISE_MOD7[x])
    
#endif 

#endif // LORE_LUTS_H