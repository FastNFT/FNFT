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



#ifdef FNFT_ENABLE_SHORT_NAMES
#define manakov_discretization_degree(...) fnft__manakov_discretization_degree(__VA_ARGS__)
#endif

#endif