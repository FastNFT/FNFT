/*
* This file is part of FNFT.  
*                                                                  
* FNFT is free software; you can redistribute it and/or
* modify it under the terms of the version 2 of the GNU General
* Public License as published by the Free Software Foundation.
*
* FNFT is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*                                                                      
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*
* Contributors:
* Lianne de Vries (TU Delft) 2021
*/
#ifndef FNFT__MANAKOV_DISCRETIZATION_H
#define FNFT__MANAKOV_DISCRETIZATION_H

#include "fnft__akns_discretization_t.h"
#include "fnft__errwarn.h"
#include "fnft__misc.h"



FNFT_UINT fnft__manakov_discretization_degree(fnft__akns_discretization_t
        discretization);

FNFT_UINT fnft__manakov_discretization_upsampling_factor(fnft__akns_discretization_t discretization);

FNFT_UINT fnft__manakov_discretization_method_order(fnft__akns_discretization_t discretization);

FNFT_INT fnft__manakov_discretization_lambda_to_z(const FNFT_UINT n, const FNFT_REAL eps_t,
	FNFT_COMPLEX* const vals, fnft__akns_discretization_t discretization);

FNFT_INT fnft__manakov_discretization_phase_factor_rho(const FNFT_REAL eps_t, const FNFT_REAL T1,
	FNFT_REAL* const phase_factor_rho, fnft__akns_discretization_t discretization);

FNFT_INT fnft__manakov_discretization_phase_factor_a(const FNFT_REAL eps_t, const FNFT_UINT D, FNFT_REAL const* const T,
	FNFT_REAL* const phase_factor_a, fnft__akns_discretization_t discretization);

FNFT_INT fnft__manakov_discretization_phase_factor_b(const FNFT_REAL eps_t, const FNFT_UINT D, FNFT_REAL const* const T,
	FNFT_REAL* const phase_factor_b, fnft__akns_discretization_t discretization);

FNFT_REAL fnft__manakov_discretization_boundary_coeff(fnft__akns_discretization_t discretization);

FNFT_UINT fnft__manakov_discretization_preprocess_signal(const FNFT_UINT D, FNFT_COMPLEX const* const q1,
		FNFT_COMPLEX const* const q2, FNFT_REAL const eps_t, const FNFT_INT kappa,
	FNFT_UINT* const Dsub_ptr, FNFT_COMPLEX** q1_preprocessed_ptr, FNFT_COMPLEX** q2_preprocessed_ptr,
	FNFT_UINT* const first_last_index, akns_discretization_t discretization);



#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_discretization_degree(...) fnft__manakov_discretization_degree(__VA_ARGS__)
#define manakov_discretization_upsampling_factor(...) fnft__manakov_discretization_upsampling_factor(__VA_ARGS__)
#define manakov_discretization_method_order(...) fnft__manakov_discretization_method_order(__VA_ARGS__)
#define manakov_discretization_lambda_to_z(...) fnft__manakov_discretization_lambda_to_z(__VA_ARGS__)
#define manakov_discretization_phase_factor_rho(...) fnft__manakov_discretization_phase_factor_rho(__VA_ARGS__)
#define manakov_discretization_phase_factor_a(...) fnft__manakov_discretization_phase_factor_a(__VA_ARGS__)
#define manakov_discretization_phase_factor_b(...) fnft__manakov_discretization_phase_factor_b(__VA_ARGS__)
#define manakov_discretization_boundary_coeff(...) fnft__manakov_discretization_boundary_coeff(__VA_ARGS__)
#define manakov_discretization_preprocess_signal(...) fnft__manakov_discretization_preprocess_signal(__VA_ARGS__)
#endif

#endif