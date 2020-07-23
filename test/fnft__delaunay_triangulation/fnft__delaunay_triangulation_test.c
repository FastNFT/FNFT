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
* Shrinivas Chimmalgi (TU Delft) 2020.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include <stdio.h>

#include "fnft__delaunay_triangulation.h"

static INT delaunay_triangulation_test()
{
    INT ret_code = SUCCESS;
    UINT NrOfTriangles = 42, i;
    UINT NodeOfTriangles1[42] = {4, 5, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 12, 
    13, 13, 14, 14, 15, 16, 16, 17, 17, 18, 18, 20, 21, 21, 22, 22, 23, 24,
    24, 25, 25, 26, 26, 8, 4, 16, 12, 24, 20};
    UINT NodeOfTriangles2[42] = {8, 9, 9, 10, 10, 11, 13, 13, 14, 14, 15, 
    15, 16, 17, 17, 18, 18, 19, 21, 21, 22, 22, 23, 23, 24, 25, 25, 26, 26, 
    27, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34};
    UINT NodeOfTriangles3[42] = {5, 8, 6, 9, 7, 10, 12, 9, 13, 10, 14, 11, 
    13, 16, 14, 17, 15, 18, 20, 17, 21, 18, 22, 19, 21, 24, 22, 25, 23, 26, 
    28, 25, 29, 26, 30, 27, 12, 8, 20, 16, 28, 24};
    
    REAL NewNodesCoordX[31] = {
      -5.000000000000000e+00,
	  -5.000000000000000e+00,
	  -5.000000000000000e+00,
	  -5.000000000000000e+00,
	  -3.333333333333333e+00,
	  -3.333333333333333e+00,
	  -3.333333333333333e+00,
	  -3.333333333333333e+00,
	  -1.666666666666667e+00,
	  -1.666666666666667e+00,
	  -1.666666666666667e+00,
	  -1.666666666666667e+00,
	  0.000000000000000e+00,
	  0.000000000000000e+00,
	  0.000000000000000e+00,
	  0.000000000000000e+00,
	  1.666666666666667e+00,
	  1.666666666666667e+00,
	  1.666666666666667e+00,
	  1.666666666666667e+00,
	  3.333333333333334e+00,
	  3.333333333333334e+00,
	  3.333333333333334e+00,
	  3.333333333333334e+00,
	  5.000000000000000e+00,
	  5.000000000000000e+00,
	  5.000000000000000e+00,
	  5.000000000000000e+00,
	  -3.333333333333333e+00,
	  0.000000000000000e+00,
	  3.333333333333334e+00
    };

    REAL NewNodesCoordY[31] = {
	  0.000000000000000e+00,
	  1.666666666666667e+00,
	  3.333333333333333e+00,
	  5.000000000000000e+00,
	  8.333333333333334e-01,
	  2.500000000000000e+00,
	  4.166666666666667e+00,
	  5.000000000000000e+00,
	  0.000000000000000e+00,
	  1.666666666666667e+00,
	  3.333333333333333e+00,
	  5.000000000000000e+00,
	  8.333333333333334e-01,
	  2.500000000000000e+00,
	  4.166666666666667e+00,
	  5.000000000000000e+00,
	  0.000000000000000e+00,
	  1.666666666666667e+00,
	  3.333333333333333e+00,
	  5.000000000000000e+00,
	  8.333333333333334e-01,
	  2.500000000000000e+00,
	  4.166666666666667e+00,
	  5.000000000000000e+00,
	  0.000000000000000e+00,
	  1.666666666666667e+00,
	  3.333333333333333e+00,
	  5.000000000000000e+00,
	  0.000000000000000e+00,
	  0.000000000000000e+00,
	  0.000000000000000e+00
    };
    UINT NrOfNewNodesCoord = 31;
    
    fnft__delaunay_triangulation_data_t * DT_ptr = NULL;
    fnft__delaunay_triangulation_data_t DT;
    DT_ptr = &DT;
    DT_ptr->NrOfNodes = 0;
    DT_ptr->NodesMax = 0;
    DT_ptr->NodeOfTriangles1 = NULL;
    DT_ptr->NodeOfTriangles2 = NULL;
    DT_ptr->NodeOfTriangles3 = NULL;
    DT_ptr->EdgesOfTriangles1 = NULL;
    DT_ptr->EdgesOfTriangles2 = NULL;
    DT_ptr->EdgesOfTriangles3 = NULL;
    DT_ptr->Edges1 = NULL;
    DT_ptr->Edges2 = NULL;
    DT_ptr->NodesCoordX = NULL;
    DT_ptr->NodesCoordY = NULL;
    DT_ptr->XCoordOfCCOfTriangles = NULL;
    DT_ptr->YCoordOfCCOfTriangles = NULL;
    DT_ptr->RadiusOfCCOfTriangles = NULL;
    DT_ptr->StatusOfTriangles = NULL;
    DT_ptr->CheckEdge = NULL;
    UINT status_flag = 0;
    
    ret_code = delaunay_triangulation(DT_ptr, 100, NrOfNewNodesCoord,
            NewNodesCoordX, NewNodesCoordY, &status_flag);
    if (ret_code != SUCCESS)
        goto leave_fun;
    
    if ((status_flag != 0) ||
            (DT_ptr->NrOfTriangles != 42))
        return E_TEST_FAILED;
    
#ifdef DEBUG
    printf("NrOfNodes=%ld \n",DT_ptr->NrOfNodes);
    printf("NodesMax=%ld \n",DT_ptr->NodesMax);
    printf("NrOfTriangles=%ld \n",DT_ptr->NrOfTriangles);
    printf("NrOfEdges=%ld \n",DT_ptr->NrOfEdges);
    for (i=0; i<NrOfTriangles; i++)
        printf("[%ld,%ld,%ld] \n",DT_ptr->NodeOfTriangles1[i],DT_ptr->NodeOfTriangles2[i],DT_ptr->NodeOfTriangles3[i]);
#endif
    for (i=0; i<NrOfTriangles; i++){
        if ((NodeOfTriangles1[i] != DT_ptr->NodeOfTriangles1[i]) ||
                (NodeOfTriangles2[i] != DT_ptr->NodeOfTriangles2[i]) ||
                (NodeOfTriangles3[i] != DT_ptr->NodeOfTriangles3[i]))
            return E_TEST_FAILED;
    }
    
    leave_fun:
        free(DT_ptr->NodeOfTriangles1);
        free(DT_ptr->NodeOfTriangles2);
        free(DT_ptr->NodeOfTriangles3);
        free(DT_ptr->EdgesOfTriangles1);
        free(DT_ptr->EdgesOfTriangles2);
        free(DT_ptr->EdgesOfTriangles3);
        free(DT_ptr->Edges1);
        free(DT_ptr->Edges2);
        free(DT_ptr->NodesCoordX);
        free(DT_ptr->NodesCoordY);
        free(DT_ptr->XCoordOfCCOfTriangles);
        free(DT_ptr->YCoordOfCCOfTriangles);
        free(DT_ptr->RadiusOfCCOfTriangles);
        free(DT_ptr->StatusOfTriangles);
        free(DT_ptr->CheckEdge);
        return ret_code;
}

INT main()
{
    if ( delaunay_triangulation_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
