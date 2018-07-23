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
* Sander Wahls (TU Delft) 2018.
*/

#include "mex.h"
#include "matrix.h"
#include "fnft.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs > 0)
        mexErrMsgTxt("No inputs expected.");

    if (nlhs == 0) {
        mexPrintf("%d.%d.%d%s\n", FNFT_VERSION_MAJOR, FNFT_VERSION_MINOR,
                  FNFT_VERSION_PATCH, FNFT_VERSION_SUFFIX);
        return;
    }

    if (nlhs != 4)
        mexErrMsgTxt("Zero or four outputs expected.");

    plhs[0] = mxCreateDoubleScalar(FNFT_VERSION_MAJOR);
    plhs[1] = mxCreateDoubleScalar(FNFT_VERSION_MINOR);
    plhs[2] = mxCreateDoubleScalar(FNFT_VERSION_PATCH);
    plhs[3] = mxCreateString(FNFT_VERSION_SUFFIX);

    return;
}
