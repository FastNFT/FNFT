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
 * Sander Wahls (TU Delft) 2020.
 */

#define FNFT_ENABLE_SHORT_NAMES
#include "fnft__nse_finvscatter.h"

#define FNFT__AKNS_FINVSCATTER_Q_FROM_R(R) (-kappa*CONJ(R))
#define FNFT__AKNS_FINVSCATTER_Q_FROM_R_PREFIX(s) fnft_nse_##s
#define FNFT__AKNS_FINVSCATTER_Q_FROM_R__PREFIX(s) fnft__nse_##s
#define FNFT__AKNS_FINVSCATTER_Q_FROM_R_HAS_KAPPA
#define FNFT__AKNS_FINVSCATTER_Q_FROM_R_HAS_MODAL

#include "fnft__akns_finvscatter_q_from_r.inc"

#undef FNFT__AKNS_FINVSCATTER_Q_FROM_R
#undef FNFT__AKNS_FINVSCATTER_Q_FROM_R_PREFIX
#undef FNFT__AKNS_FINVSCATTER_Q_FROM_R__PREFIX
#undef FNFT__AKNS_FINVSCATTER_Q_FROM_R_HAS_KAPPA
#undef FNFT__AKNS_FINVSCATTER_Q_FROM_R_HAS_MODAL
