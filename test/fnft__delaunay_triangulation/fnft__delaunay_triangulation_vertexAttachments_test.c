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
    UINT NrOfNodes = 5;  
    UINT Nodes[5] = {6, 10, 17, 4, 8};
    UINT ArrayOfCandidateElements_exact[23] = {2, 3, 4, 3, 4, 5, 9, 10, 11, 
    13, 14 , 15 , 19, 20, 21, 0, 37, 0, 1, 6, 7, 36, 37};
    UINT NrOfCandidateElements_exact = 23;
    UINT *ArrayOfCandidateElements = NULL;
    ArrayOfCandidateElements = malloc(3*NrOfTriangles*sizeof(UINT));
    if (ArrayOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    UINT NrOfCandidateElements = 0;
          
    ret_code = delaunay_triangulation_vertexAttachments(NrOfTriangles,
        NodeOfTriangles1, NodeOfTriangles2, NodeOfTriangles3,
        NrOfNodes, Nodes, &NrOfCandidateElements, ArrayOfCandidateElements);   
        
    if (ret_code != SUCCESS)
        goto leave_fun;
    
    if ((ArrayOfCandidateElements == NULL) ||
            (NrOfCandidateElements != NrOfCandidateElements_exact))
        return E_TEST_FAILED;
    
#ifdef DEBUG
    printf("NrOfCandidateElements=%ld \n",NrOfCandidateElements);
    for (i=0; i<NrOfCandidateElements; i++)
        printf("%ld \n",ArrayOfCandidateElements[i]);
#endif
    for (i=0; i<NrOfCandidateElements_exact; i++){
        if (ArrayOfCandidateElements[i] != ArrayOfCandidateElements_exact[i])
            return E_TEST_FAILED;
    }
    
    leave_fun:
        free(ArrayOfCandidateElements);
        return ret_code;
}

INT main()
{
    if ( delaunay_triangulation_test() != SUCCESS )
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
