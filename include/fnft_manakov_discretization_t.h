/*
This file specifies the types for the manakov equation
*/

#ifndef FNFT__MANAKOV_DISCRETIZATION_T_H    //header guard
#define FNFT__MANAKOV_DISCRETIZATION_T_H

#include "fnft.h"

typedef enum {
    fnft_manakov_discretization_2SPLIT2_MODAL,
    fnft_manakov_discretization_BO,
    fnft_manakov_discretization_2SPLIT1A,
    fnft_manakov_discretization_2SPLIT1B,
    fnft_manakov_discretization_2SPLIT2A,
    fnft_manakov_discretization_2SPLIT2B,
    fnft_manakov_discretization_2SPLIT2S,
    fnft_manakov_discretization_2SPLIT3A,
    fnft_manakov_discretization_2SPLIT3B,
    fnft_manakov_discretization_2SPLIT3S,
    fnft_manakov_discretization_2SPLIT4A,
    fnft_manakov_discretization_2SPLIT4B,
    fnft_manakov_discretization_2SPLIT5A,
    fnft_manakov_discretization_2SPLIT5B,
    fnft_manakov_discretization_2SPLIT6A,
    fnft_manakov_discretization_2SPLIT6B,
    fnft_manakov_discretization_2SPLIT7A,
    fnft_manakov_discretization_2SPLIT7B,
    fnft_manakov_discretization_2SPLIT8A,
    fnft_manakov_discretization_2SPLIT8B,
    fnft_manakov_discretization_4SPLIT4A,
    fnft_manakov_discretization_4SPLIT4B,
    fnft_manakov_discretization_CF4_2,
    fnft_manakov_discretization_CF4_3,
    fnft_manakov_discretization_CF5_3,
    fnft_manakov_discretization_CF6_4,
    fnft_manakov_discretization_ES4,
    fnft_manakov_discretization_TES4
} fnft_manakov_discretization_t;

// TODO: add short names code, see nse file

#endif
