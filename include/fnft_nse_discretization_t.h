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
* Sander Wahls (TU Delft) 2017.
*/

/**
 * @file fnft_nse_discretization_t.h
 * @brief Lists discretizations for the nonlinear Schroedinger equation.
 * @ingroup fnft
 */
#ifndef FNFT_NSE_DISCRETIZATION_T_H
#define FNFT_NSE_DISCRETIZATION_T_H

#include "fnft.h"

/**
 * Enum that specifies discretizations used to compute nonlinear Fourier
 * transforms for the nonlinear Schroedinger equation. Used in
 * \link fnft_nsev_opts_t \endlink.\n \n
 * fnft_nse_discretization_2SPLIT2A_MODAL: It is the normalized Ablowitz-Ladik discretization from 
 * Wahls and Poor,<a href="http://dx.doi.org/10.1109/ICASSP.2013.6638772">&quot;Introducing the fast nonlinear Fourier transform,&quot;</a> Proc. ICASSP 2013.\n \n
 * fnft_nse_discretization_2SPLIT2A : It has been taken from \n \n
 * fnft_nse_discretization_2SPLIT4A and fnft_nse_discretization_2SPLIT4B : They are equivalent forms of
 * Eq. 20 in Prins and Wahls, &quot;Higher order exponential splittings for the fast non-linear Fourier transform of the KdV equation,&quot; to appear in Proc. ICASSP 2018. \n \n
 * fnft_nse_discretization_BO : It has been taken from Boffetta and Osborne, <a href="https://doi.org/10.1016/0021-9991(92)90370-E">&quot;Computation of the direct scattering transform for the nonlinear Schroedinger  equation,&quot;</a> J. Comput. Phys. 102(2), 1992. It is supported by \link fnft__nse_scatter.h \endlink.\n 
 * Other discretizations are supported by \link fnft__nse_fscatter.h \endlink.\n
 * @ingroup data_types
 */
typedef enum {
    fnft_nse_discretization_2SPLIT2_MODAL,
    fnft_nse_discretization_BO,
    fnft_nse_discretization_2SPLIT1A,
    fnft_nse_discretization_2SPLIT1B,
    fnft_nse_discretization_2SPLIT2A,
    fnft_nse_discretization_2SPLIT2B,
    fnft_nse_discretization_2SPLIT2S,
    fnft_nse_discretization_2SPLIT3A,
    fnft_nse_discretization_2SPLIT3B,
    fnft_nse_discretization_2SPLIT3S,
    fnft_nse_discretization_2SPLIT4A,
    fnft_nse_discretization_2SPLIT4B,
    fnft_nse_discretization_2SPLIT5A,
    fnft_nse_discretization_2SPLIT5B,
    fnft_nse_discretization_2SPLIT6A,
    fnft_nse_discretization_2SPLIT6B,
    fnft_nse_discretization_2SPLIT7A,
    fnft_nse_discretization_2SPLIT7B,
    fnft_nse_discretization_2SPLIT8A,
    fnft_nse_discretization_2SPLIT8B
} fnft_nse_discretization_t;

#ifdef FNFT_ENABLE_SHORT_NAMES
#define nse_discretization_2SPLIT2_MODAL fnft_nse_discretization_2SPLIT2_MODAL
#define nse_discretization_BO fnft_nse_discretization_BO
#define nse_discretization_t fnft_nse_discretization_t
#define nse_discretization_2SPLIT1A fnft_nse_discretization_2SPLIT1A
#define nse_discretization_2SPLIT1B fnft_nse_discretization_2SPLIT1B
#define nse_discretization_2SPLIT2A fnft_nse_discretization_2SPLIT2A
#define nse_discretization_2SPLIT2B fnft_nse_discretization_2SPLIT2B
#define nse_discretization_2SPLIT2S fnft_nse_discretization_2SPLIT2S
#define nse_discretization_2SPLIT3A fnft_nse_discretization_2SPLIT3A
#define nse_discretization_2SPLIT3B fnft_nse_discretization_2SPLIT3B
#define nse_discretization_2SPLIT3S fnft_nse_discretization_2SPLIT3S
#define nse_discretization_2SPLIT4A fnft_nse_discretization_2SPLIT4A
#define nse_discretization_2SPLIT4B fnft_nse_discretization_2SPLIT4B
#define nse_discretization_2SPLIT5A fnft_nse_discretization_2SPLIT5A
#define nse_discretization_2SPLIT5B fnft_nse_discretization_2SPLIT5B
#define nse_discretization_2SPLIT6A fnft_nse_discretization_2SPLIT6A
#define nse_discretization_2SPLIT6B fnft_nse_discretization_2SPLIT6B
#define nse_discretization_2SPLIT7A fnft_nse_discretization_2SPLIT7A
#define nse_discretization_2SPLIT7B fnft_nse_discretization_2SPLIT7B
#define nse_discretization_2SPLIT8A fnft_nse_discretization_2SPLIT8A
#define nse_discretization_2SPLIT8B fnft_nse_discretization_2SPLIT8B
#endif

#endif
