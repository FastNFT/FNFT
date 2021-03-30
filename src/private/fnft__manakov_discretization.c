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
* Contributor:
* Lianne de Vries (TU Delft) 2021.
*/


#define FNFT_ENABLE_SHORT_NAMES

#include "fnft__manakov_discretization.h"

/**
 * This routine returns the max degree d of the polynomials in a single
 * scattering matrix or zero if the discretization is unknown.
 * 
 */
UINT fnft__manakov_discretization_degree(akns_discretization_t
        discretization)
{
    switch (discretization) {
        case akns_discretization_2SPLIT3A:
        case akns_discretization_2SPLIT3B:
            return 6;
        case akns_discretization_4SPLIT4A:
        case akns_discretization_4SPLIT4B:
        case akns_discretization_FTES4_4A:
        case akns_discretization_FTES4_4B:
            return 8;

                    
        default: // Unknown discretization
            return 0;
    }
}



