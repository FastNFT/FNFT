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