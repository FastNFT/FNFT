#ifndef FNFT_MANAKOV_H
#define FNFT_MANAKOV_H

#include "fnft__manakov_fscatter.h"
// #include "fnft__nse_scatter.h"
#include "fnft__misc.h" // for l2norm
#include <string.h> // for memcpy
#include <stdio.h>
#include "fnft__errwarn.h"
#include "fnft__poly_roots_fasteigen.h"
#include "fnft__poly_chirpz.h"
#include "fnft_manakov_discretization_t.h"


typedef enum{
    fnft_manakov_bsfilt_NONE,
    fnft_manakov_bsfilt_BASIC,
    fnft_manakov_bsfilt_FULL
} fnft_manakov_bsfilt_t;


typedef enum{
    fnft_manakov_bsloc_FAST_EIGENVALUE,
    fnft_manakov_bsloc_NEWTON,
    fnft_manakov_bsloc_SUBSAMPLE_AND_REFINE
} fnft_manakov_bsloc_t;

typedef enum{
    fnft_manakov_dstype_NORMING_CONSTANTS,
    fnft_manakov_dstype_RESIDUES,
    fnft_manakov_dstype_BOTH
} fnft_manakov_dstype_t;

typedef enum {
	fnft_manakov_cstype_REFLECTION_COEFFICIENT,
	fnft_manakov_cstype_AB,
	fnft_manakov_cstype_BOTH
} fnft_manakov_cstype_t;

typedef struct {
	fnft_manakov_bsfilt_t bound_state_filtering;
	fnft_manakov_bsloc_t bound_state_localization;
	FNFT_UINT niter;
	FNFT_UINT Dsub;
	fnft_manakov_dstype_t discspec_type;
	fnft_manakov_cstype_t contspec_type;
	FNFT_INT normalization_flag;
	fnft_manakov_discretization_t discretization;
	FNFT_UINT richardson_extrapolation_flag;
} fnft_manakov_opts_t;

fnft_manakov_opts_t fnft_manakov_default_opts();

FNFT_INT fnft_manakov(const FNFT_UINT D, FNFT_COMPLEX const* const q1, FNFT_COMPLEX const* const q2,
	FNFT_REAL const* const T, const FNFT_UINT M,
	FNFT_COMPLEX* const contspec, FNFT_REAL const* const XI,
	FNFT_UINT* const K_ptr, FNFT_COMPLEX* const bound_states,
	FNFT_COMPLEX* const normconsts_or_residues, const FNFT_INT kappa,
	fnft_manakov_opts_t* opts);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_bsfilt_NONE fnft_manakov_bsfilt_NONE
#define manakov_bsfilt_BASIC fnft_manakov_bsfilt_BASIC
#define manakov_bsfilt_FULL fnft_manakov_bsfilt_FULL
#define manakov_bsloc_FAST_EIGENVALUE fnft_manakov_bsloc_FAST_EIGENVALUE
#define manakov_bsloc_NEWTON fnft_manakov_bsloc_NEWTON
#define manakov_bsloc_SUBSAMPLE_AND_REFINE fnft_manakov_bsloc_SUBSAMPLE_AND_REFINE
#define manakov_dstype_NORMING_CONSTANTS fnft_manakov_dstype_NORMING_CONSTANTS
#define manakov_dstype_RESIDUES fnft_manakov_dstype_RESIDUES
#define manakov_dstype_BOTH fnft_manakov_dstype_BOTH
#define manakov_cstype_REFLECTION_COEFFICIENT fnft_manakov_cstype_REFLECTION_COEFFICIENT
#define manakov_cstype_AB fnft_manakov_cstype_AB
#define manakov_cstype_BOTH fnft_manakov_cstype_BOTH
#endif

#endif      // FNFT_MANAKOV_H