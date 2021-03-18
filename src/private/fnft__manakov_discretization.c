// Not sure yet what functions should be in this file and if it is necessary,
// but because the methods for 2x2 all use a change of base trick to reduce the degree,
// we need a different function to determine the degree for the 3x3 system methods


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

                    
        default: // Unknown discretization
            return 0;
    }
}



