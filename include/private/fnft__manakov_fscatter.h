

#ifndef FNFT__MANAKOV_FSCATTER_H
#define FNFT__MANAKOV_FSCATTER_H

#include "fnft__akns_discretization.h"
#include "fnft__manakov_discretization.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"


FNFT_UINT fnft__manakov_fscatter_numel(FNFT_UINT D,
	fnft__akns_discretization_t discretization);

FNFT_INT fnft__manakov_fscatter(const UINT D, COMPLEX const* const q, const UINT kappa,
	REAL eps_t, COMPLEX* const result, UINT* const deg_ptr,
	INT* const W_ptr, akns_discretization_t const discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_fscatter_numel(...) fnft__manakov_fscatter_numel(__VA_ARGS__)
#define manakov_fscatter(...) fnft__manakov_fscatter(__VA_ARGS__)
# endif


#endif // !FNFT__MANAKOV_FSCATTER_H



                                    