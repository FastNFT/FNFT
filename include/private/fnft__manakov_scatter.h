
#ifndef FNFT__MANAKOV_SCATTER_H
#define FNFT__MANAKOV_SCATTER_H

#include "fnft__manakov_discretization.h"
//#include "fnft_manakov_discretization_t.h"
//#include "fnft_numtypes.h"
#include <stdio.h>
#include <string.h> // for memcpy
#include "fnft__errwarn.h"
#include "fnft__misc.h"// for square_matrix_mult


FNFT_INT fnft__manakov_scatter_matrix(FNFT_UINT const D,
                        FNFT_COMPLEX const * const q1,
                        FNFT_COMPLEX const * const q2,
                        FNFT_REAL const eps_t,
                        FNFT_UINT const K,
                        FNFT_COMPLEX const * const lambda,
                        FNFT_INT const kappa,
                        FNFT_COMPLEX * const result,
                        fnft_manakov_discretization_t const discretization);


#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_scatter_matrix(...) fnft__manakov_scatter_matrix(__VA_ARGS__)
#endif

#endif
