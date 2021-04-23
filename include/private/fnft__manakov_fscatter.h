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

#ifndef FNFT__MANAKOV_FSCATTER_H
#define FNFT__MANAKOV_FSCATTER_H

#include "fnft__manakov_discretization.h"
#include "fnft__misc.h"
#include "fnft__errwarn.h"
#include "fnft__poly_fmult.h"


FNFT_UINT fnft__manakov_fscatter_numel(FNFT_UINT D,
	fnft_manakov_discretization_t discretization);

FNFT_INT fnft__manakov_fscatter(const UINT D, COMPLEX const* const q1, COMPLEX const* const q2, const UINT kappa,
	REAL eps_t, COMPLEX* const result, UINT* const deg_ptr,
	INT* const W_ptr, manakov_discretization_t const discretization);

#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_fscatter_numel(...) fnft__manakov_fscatter_numel(__VA_ARGS__)
#define manakov_fscatter(...) fnft__manakov_fscatter(__VA_ARGS__)
# endif


#endif // !FNFT__MANAKOV_FSCATTER_H



                                    