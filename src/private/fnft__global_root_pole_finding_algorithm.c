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
* Sander Wahls (TU Delft) 2022.
*
* This code is based on a Matlab implementation of the algorithm by
*
* Piotr Kowalczyk (Gdansk University of Technology), 2018-2019
*
* that is available at https://github.com/PioKow/GRPF under the terms
* of the MIT license reproduced below:
*
* Copyright (c) 2018 Gdansk University of Technology
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*/

#define FNFT_ENABLE_SHORT_NAMES

#include "fnft.h"
#include "fnft__poly_eval.h"
#include "fnft__errwarn.h"
#include "fnft__delaunay_triangulation.h"

static UINT FindNextNode(REAL *NodesCoordX, REAL *NodesCoordY, UINT PrevNode, 
        UINT RefNode, UINT *TempNodes, UINT NrOfTempNodes)
{
// %  FindNextNode: finds the next node in the candidate region boudary
// %                process. The next one (after the reference one) is picked
// %                from the fixed set of nodes.
// %
// % INPUTS
// %
// %  NodesCoord     : nodes coordinates
// %  PrevNode       : previous node
// %  RefNode        : reference (current) node
// %  TempNodes      : set of nodes
// %
// % OUTPUTS
// %
// %  Index          : index of the next node
// %

    UINT Index = 0, j;
    REAL min_Phi = 10;
    REAL PrevNodeCoordX = NodesCoordX[PrevNode];
    REAL PrevNodeCoordY = NodesCoordY[PrevNode];
    REAL RefNodeCoordX = NodesCoordX[RefNode];
    REAL RefNodeCoordY = NodesCoordY[RefNode];
    REAL PRX, PRY, TRX, TRY, LenPR, LenTR, DotProd, Phi;
    for (j=0; j<NrOfTempNodes; j++){
        PRX = PrevNodeCoordX-RefNodeCoordX;
        PRY = PrevNodeCoordY-RefNodeCoordY;
        TRX = NodesCoordX[TempNodes[j]] - RefNodeCoordX;
        TRY = NodesCoordY[TempNodes[j]] - RefNodeCoordY;
        LenPR = SQRT(PRX*PRX + PRY*PRY);
        LenTR = SQRT(TRX*TRX + TRY*TRY);
        DotProd = PRX*TRX + PRY*TRY;
        Phi = ACOS(DotProd/(LenPR*LenTR));
        if ((PRX*TRY-PRY*TRX)<0)
            Phi = 2*PI - Phi;
        if (Phi < min_Phi){
            min_Phi = Phi;
            Index = j;
        }
    }
    
    return Index;
}

INT fnft__global_root_pole_finding_algorithm(
        UINT * const K_roots_ptr,
        UINT * roots_multiplicity,
        COMPLEX * roots,
        UINT * const K_poles_ptr,
        UINT * poles_multiplicity,
        COMPLEX * poles,
        INT fun (UINT, COMPLEX *, COMPLEX *, void *),
        void * params_ptr,
        UINT NodesMax,
        const REAL bounding_box_local[4],
        const REAL Tol,
        const UINT niter)
{
    INT ret_code = SUCCESS;
    UINT i, j;
    
    UINT NodesShift = 0;
    UINT NrOfElements = 0;
    UINT NrOfValidEdges = 0;
    UINT NrOfCandidateEdges = 0;
    UINT NrOfCandidateNodes = 0;
    UINT NrOfCandidateElements = 0;
    UINT *ArrayOfCandidateElements = NULL;
    UINT NrOfFirstZoneCandidateElements = 0;
    UINT NrOfSecondZoneCandidateElements = 0;
    REAL XCoordOfTempEdgeNode1, XCoordOfTempEdgeNode2;
    REAL YCoordOfTempEdgeNode1, YCoordOfTempEdgeNode2;
    REAL XTempNodeCoord, YTempNodeCoord, TempEdgeLength, DistNodes;
    UINT AddCoordFlag;
    UINT *IDOfCandidateElements = NULL;
    UINT *MultiplicationOfTempEdges = NULL;
    UINT *ContourEdges1 = NULL;
    UINT *ContourEdges2 = NULL;
    UINT *NrOfNodesOfRegions = NULL, *Regions = NULL, *IndexOfNextEdge = NULL, *TempNodes = NULL, *SideFlag = NULL;
    INT *z_m = NULL;
    COMPLEX *z = NULL;
    REAL *NewNodesCoordX = NULL;
    REAL *NewNodesCoordY = NULL;
    UINT NrOfNewNodes = 0;
    UINT NrOfNodes = 0;
    REAL *NodesCoordX = NULL;
    REAL *NodesCoordY = NULL;
    COMPLEX * FuntionValues = NULL;
    COMPLEX * ComplexNodes = NULL;
    UINT * Quadrants = NULL;
    REAL * Angles = NULL;
    UINT *ValidEdges = NULL;
    UINT * Elements1 = NULL, * Elements2 = NULL, * Elements3 = NULL;
    UINT *ValidEdges1 = NULL, *ValidEdges2 = NULL;
    UINT *CandidateEdges1 = NULL, *CandidateEdges2 = NULL;
    UINT *ElementsCount = NULL, *CandidateNodes = NULL, * Temp = NULL;
    UINT *FlagOfFirstZoneCandidateElements = NULL;
    UINT *FlagOfSecondZoneCandidateElements = NULL;
    UINT *TempExtraEdges1 = NULL, *TempExtraEdges2 = NULL;
    UINT DT_built = 0;

    COMPLEX * ComplexPoints = NULL;
    UINT * NewQuadrants = NULL;
    REAL * NewAngles = NULL;
    INT * contour_orientation = NULL;
    
    
    UINT * NodeOfTriangles1 = NULL;
    UINT * NodeOfTriangles2 = NULL;
    UINT * NodeOfTriangles3 = NULL;
    UINT * EdgesOfTriangles1 = NULL;
    UINT * EdgesOfTriangles2 = NULL;
    UINT * EdgesOfTriangles3 = NULL;
    UINT * Edges1 = NULL;
    UINT * Edges2 = NULL;
    UINT NrOfEdges = 0;
    UINT NrOfTriangles = 0;
    UINT NrOfNodesToRemove = 0;

    ret_code =  delaunay_triangulation_rect_dom(bounding_box_local[0], bounding_box_local[1],
            bounding_box_local[2], bounding_box_local[3],
            (bounding_box_local[3] - bounding_box_local[2])/5,
            &NrOfNewNodes, &NewNodesCoordX, &NewNodesCoordY);
    CHECK_RETCODE(ret_code, leave_fun);
//     printf("NrOfNewNodes=%ld \n",NrOfNewNodes);

    fnft__delaunay_triangulation_data_t * DT_ptr = NULL;
    fnft__delaunay_triangulation_data_t DT;
    DT_ptr = &DT;
    DT_ptr->NrOfNodes = 0;
    DT_ptr->NodesMax = 0;    
    UINT status_flag = 0;
    // Setting all poniters to NULL so that they can be
    // freed if allocation fails for some pointer
    // before it.
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
    DT_built = 1;
    
    NodesCoordX = malloc(NodesMax * sizeof(REAL));
    NodesCoordY = malloc(NodesMax * sizeof(REAL));
    if (NodesCoordX == NULL || NodesCoordX == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    FuntionValues = malloc(NodesMax * sizeof(COMPLEX));
    ComplexNodes = malloc(NodesMax * sizeof(COMPLEX));
    Quadrants = malloc(NodesMax * sizeof(UINT));
    Angles = malloc(NodesMax * sizeof(REAL));
    if (FuntionValues == NULL || Quadrants == NULL 
            || ComplexNodes == NULL || Angles == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ValidEdges = malloc(1*sizeof(UINT));
    if (ValidEdges == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    Elements1 = malloc(1*sizeof(UINT));
    Elements2 = malloc(1*sizeof(UINT));
    Elements3 = malloc(1*sizeof(UINT));
    if (Elements1 == NULL || Elements2 == NULL|| Elements3 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ValidEdges1 = malloc(1*sizeof(UINT));
    ValidEdges2 = malloc(1*sizeof(UINT));
    if (ValidEdges1 == NULL || ValidEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ArrayOfCandidateElements = malloc(1*sizeof(UINT));        
    CandidateEdges1 = malloc(1*sizeof(UINT));
    CandidateEdges2 = malloc(1*sizeof(UINT));
    if (CandidateEdges1 == NULL || CandidateEdges2 == NULL
            || ArrayOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }    
    
    Temp = malloc(NodesMax*sizeof(UINT));
     if (Temp == NULL ){
        ret_code = E_NOMEM;
        goto leave_fun;
    }    
   
    CandidateNodes = malloc(1*sizeof(UINT));
    if (CandidateNodes == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    ElementsCount = malloc(1*sizeof(UINT));
    if (ElementsCount == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    FlagOfFirstZoneCandidateElements = malloc(1*sizeof(UINT));
    FlagOfSecondZoneCandidateElements = malloc(1*sizeof(UINT));
    if (FlagOfFirstZoneCandidateElements == NULL || FlagOfSecondZoneCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    TempExtraEdges1 = malloc(1*sizeof(UINT));
    TempExtraEdges2 = malloc(1*sizeof(UINT));
    if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
   
    INT converged_flag = 0;
    for (UINT iter=1; iter<=niter; iter++){
        
//         printf("iter=%ld\n",iter);
//         printf("NrOfNewNodes=%ld \n",NrOfNewNodes);
        if (NrOfNodes+NrOfNewNodes > NodesMax)
            break;
        
        NodesShift = NrOfNodes;
        NrOfNodes = NrOfNodes+NrOfNewNodes;
        
        for (i=0; i<NrOfNewNodes; i++){
            NodesCoordX[NodesShift+i] = NewNodesCoordX[i];
            NodesCoordY[NodesShift+i] = NewNodesCoordY[i];
            ComplexNodes[i] = NewNodesCoordX[i]+I*NewNodesCoordY[i];
        }
        
        ret_code = fun(NrOfNewNodes, ComplexNodes, FuntionValues+NodesShift, params_ptr);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        ret_code = misc_quadrant(NrOfNewNodes, FuntionValues+NodesShift,
                Quadrants+NodesShift, Angles+NodesShift);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
//         printf("[");        
//         for (j=0; j<NrOfNodes; j++)
//             printf("%ld, ", Quadrants[j]);
//         printf("]\n");

        ret_code = delaunay_triangulation(DT_ptr, NodesMax, NrOfNewNodes,
                NewNodesCoordX, NewNodesCoordY, &status_flag);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
        
        if (status_flag == 1){
            WARN("Duplicate point.");
            break;
        }
        if (status_flag == 2){
            WARN("Ran out of memory for triangles.");
            break;
        }
         if (status_flag == 3){
            WARN("Ran out of memory for Edges.");
            break;
        }
        
//         printf("NrOfNodes=%ld \n",DT_ptr->NrOfNodes);
//         printf("NodesMax=%ld \n",DT_ptr->NodesMax);
//         printf("NrOfTriangles=%ld \n",DT_ptr->NrOfTriangles);
//         printf("NrOfEdges=%ld \n",DT_ptr->NrOfEdges);
        
        NodeOfTriangles1 = DT_ptr->NodeOfTriangles1;
        NodeOfTriangles2 = DT_ptr->NodeOfTriangles2;
        NodeOfTriangles3 = DT_ptr->NodeOfTriangles3;
        Edges1 = DT_ptr->Edges1;
        Edges2 = DT_ptr->Edges2;
        NrOfEdges = DT_ptr->NrOfEdges;
        EdgesOfTriangles1 = DT_ptr->EdgesOfTriangles1;
        EdgesOfTriangles2 = DT_ptr->EdgesOfTriangles2;
        EdgesOfTriangles3 = DT_ptr->EdgesOfTriangles3;
        NrOfTriangles = DT_ptr->NrOfTriangles;
        NrOfNodesToRemove = DT_ptr->NrOfNodesToRemove;
        
        ValidEdges = realloc(ValidEdges,NrOfEdges*sizeof(UINT));
        if (ValidEdges == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfEdges; i++)
            ValidEdges[i] = 0;
        
        NrOfElements = NrOfTriangles;       
        for (i=0; i<NrOfTriangles; i++){
            ValidEdges[EdgesOfTriangles1[i]] = 1;
            ValidEdges[EdgesOfTriangles2[i]] = 1;
            ValidEdges[EdgesOfTriangles3[i]] = 1;
        }
        
        Elements1 = realloc(Elements1,NrOfElements*sizeof(UINT));
        Elements2 = realloc(Elements2,NrOfElements*sizeof(UINT));
        Elements3 = realloc(Elements3,NrOfElements*sizeof(UINT));
        if (Elements1 == NULL || Elements2 == NULL|| Elements3 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfTriangles; i++){
            Elements1[i] = NodeOfTriangles1[i]-NrOfNodesToRemove;
            Elements2[i] = NodeOfTriangles2[i]-NrOfNodesToRemove;
            Elements3[i] = NodeOfTriangles3[i]-NrOfNodesToRemove;
        }
                
//         for (j=0; j<NrOfElements; j++)
//             printf("[%ld,%ld,%ld]\n", Elements1[j],Elements2[j],Elements3[j]);
        
        NrOfValidEdges = 0;
        for (i=0; i<NrOfEdges; i++){
            if (ValidEdges[i] == 1)
                NrOfValidEdges = NrOfValidEdges+1;
        }
//         printf(" NrOfValidEdges=%ld \n", NrOfValidEdges);
        
        ValidEdges1 = realloc(ValidEdges1, NrOfValidEdges*sizeof(UINT));
        ValidEdges2 = realloc(ValidEdges2, NrOfValidEdges*sizeof(UINT));
        if (ValidEdges1 == NULL || ValidEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        j = 0;
        for (i=0; i<NrOfEdges; i++){
            if (ValidEdges[i] == 1){
                ValidEdges1[j] = Edges1[i]-NrOfNodesToRemove;
                ValidEdges2[j] = Edges2[i]-NrOfNodesToRemove;
//                 printf("[%ld,%ld]\n", ValidEdges1[j],ValidEdges2[j]);
                j = j+1;
            }
        }
        // CandidateEdges are edges across whom phase difference is 2
        CandidateEdges1 = realloc(CandidateEdges1, NrOfValidEdges*sizeof(UINT));
        CandidateEdges2 = realloc(CandidateEdges2, NrOfValidEdges*sizeof(UINT));
        if (CandidateEdges1 == NULL || CandidateEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        NrOfCandidateEdges = 0;
        for (i=0; i<NrOfValidEdges; i++){
            if (((Quadrants[ValidEdges1[i]]-Quadrants[ValidEdges2[i]])%4) == 2){
                CandidateEdges1[NrOfCandidateEdges] =  ValidEdges1[i];
                CandidateEdges2[NrOfCandidateEdges] =  ValidEdges2[i];
                NrOfCandidateEdges = NrOfCandidateEdges + 1;
            }
        }
        
        if  (NrOfCandidateEdges == 0){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
//         printf("NrOfCandidateEdges=%ld\n",NrOfCandidateEdges );
        
        for (i=0; i<(DT_ptr->NrOfNodes); i++)
            Temp[i] = 0;
        
        REAL MinCandidateEdgesLengths = INFINITY;
        REAL MaxCandidateEdgesLengths = 0;
        REAL CurrEdgeLength = 0;
        // Length of the CandidateEdges
        for (i=0; i<NrOfCandidateEdges; i++){
            CurrEdgeLength = SQRT((NodesCoordX[CandidateEdges2[i]]-NodesCoordX[CandidateEdges1[i]])*(NodesCoordX[CandidateEdges2[i]]-NodesCoordX[CandidateEdges1[i]])
            +(NodesCoordY[CandidateEdges2[i]]-NodesCoordY[CandidateEdges1[i]])*(NodesCoordY[CandidateEdges2[i]]-NodesCoordY[CandidateEdges1[i]]));
            if (CurrEdgeLength < MinCandidateEdgesLengths)
                MinCandidateEdgesLengths = CurrEdgeLength;
            if (CurrEdgeLength > MaxCandidateEdgesLengths)
                MaxCandidateEdgesLengths = CurrEdgeLength;
            // If a CandidateEdge is longer than Tol then its end points are marked
            if (CurrEdgeLength > Tol){
                Temp[CandidateEdges1[i]] = 1;// TODO may require to set everything to 0 first
                Temp[CandidateEdges2[i]] = 1;
            }
        }
        
        if (MaxCandidateEdgesLengths < Tol){
            //printf("Tolerence achieved in interation %ld \n", iter);
            converged_flag = 1;
            break;
        }
        
        // Collecting all the nodes which are the end points
        // Marking first and then collecting to preventing counting a node twice
        CandidateNodes = realloc(CandidateNodes, NrOfNodes*sizeof(UINT));
        if (CandidateNodes == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        NrOfCandidateNodes = 0;
        for (i=0; i<NrOfNodes; i++){
            if (Temp[i] == 1){
                CandidateNodes[NrOfCandidateNodes] = i;
                NrOfCandidateNodes = NrOfCandidateNodes+1;
             //printf("CandidateNode =%ld\n", i);
            }
        }
//         printf(" NrOfCandidateNodes =%ld\n", NrOfCandidateNodes  );

        // Listing all the triangles that the nodes in CandidateNodes are part
        // of
        NrOfCandidateElements = 0;
        ArrayOfCandidateElements= realloc(ArrayOfCandidateElements, 3*NrOfElements*sizeof(UINT));
        if (ArrayOfCandidateElements == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        ret_code =  delaunay_triangulation_vertexAttachments(NrOfElements,Elements1, Elements2, Elements3,
                NrOfCandidateNodes, CandidateNodes, &NrOfCandidateElements, ArrayOfCandidateElements);
        if (ret_code != SUCCESS){
            ret_code = E_SUBROUTINE(ret_code);
            CHECK_RETCODE(ret_code, leave_fun);
        }
//         printf(" NrOfCandidateElements =%ld\n", NrOfCandidateElements  );
//  for (i=0; i<NrOfCandidateElements; i++)
//      printf(" CandidateElement =%ld\n", ArrayOfCandidateElements[i]);
//
        // Calculating the number of times a particular triangle shows up in
        // ArrayOfCandidateElements
        ElementsCount = realloc(ElementsCount, NrOfElements*sizeof(UINT));
        if (ElementsCount == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfElements; i++)
            ElementsCount[i] = 0;
        for (i=0; i<NrOfCandidateElements; i++)
            ElementsCount[ArrayOfCandidateElements[i]] = ElementsCount[ArrayOfCandidateElements[i]]+1;
                
        NrOfFirstZoneCandidateElements = 0;
        NrOfSecondZoneCandidateElements = 0;
        
        FlagOfFirstZoneCandidateElements = realloc(FlagOfFirstZoneCandidateElements, NrOfElements*sizeof(UINT));
        FlagOfSecondZoneCandidateElements = realloc(FlagOfSecondZoneCandidateElements, NrOfElements*sizeof(UINT));
        if (FlagOfFirstZoneCandidateElements == NULL || FlagOfSecondZoneCandidateElements == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        for (i=0; i<NrOfElements; i++){
            FlagOfFirstZoneCandidateElements[i] = 0;
            FlagOfSecondZoneCandidateElements[i] = 0;
        }
        
        for (i=0; i<NrOfElements; i++){
            // If a triangle shows up more than once
            if (ElementsCount[i]>1){
                FlagOfFirstZoneCandidateElements[i] = 1;
                NrOfFirstZoneCandidateElements = NrOfFirstZoneCandidateElements+1;
            }
            // If a triangle shows up exactly once
            if (ElementsCount[i] == 1){
                FlagOfSecondZoneCandidateElements[i] = 1;
                NrOfSecondZoneCandidateElements = NrOfSecondZoneCandidateElements+1;
            }
        }
//         printf("NrOfFirstZoneCandidateElements =%ld\n", NrOfFirstZoneCandidateElements);
//         printf("NrOfSecondZoneCandidateElements =%ld\n", NrOfSecondZoneCandidateElements);
//         
        // TODO possibility to used Edges1 and Edges2 directly
        // Vectorizing FirstZoneCandidateElements
        TempExtraEdges1 = realloc(TempExtraEdges1, NrOfFirstZoneCandidateElements*3*sizeof(UINT));
        TempExtraEdges2 = realloc(TempExtraEdges2, NrOfFirstZoneCandidateElements*3*sizeof(UINT));
        if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        j = 0;
        for (i=0; i<NrOfElements; i++){
            if (FlagOfFirstZoneCandidateElements[i] == 1){
                TempExtraEdges1[j*3] = Elements1[i];
                TempExtraEdges2[j*3] = Elements2[i];
                
                TempExtraEdges1[j*3+1] = Elements2[i];
                TempExtraEdges2[j*3+1] = Elements3[i];
                
                TempExtraEdges1[j*3+2] = Elements3[i];
                TempExtraEdges2[j*3+2] = Elements1[i];
                j = j+1;
            }
        }
        
        // Creating newnodes at mid-points of the edges
        NewNodesCoordX = realloc(NewNodesCoordX, 3*(NrOfFirstZoneCandidateElements+NrOfSecondZoneCandidateElements)*sizeof(REAL));
        NewNodesCoordY = realloc(NewNodesCoordY, 3*(NrOfFirstZoneCandidateElements+NrOfSecondZoneCandidateElements)*sizeof(REAL));
        NewNodesCoordX[0] = (NodesCoordX[TempExtraEdges1[0]]+NodesCoordX[TempExtraEdges2[0]])/2;
        NewNodesCoordY[0] = (NodesCoordY[TempExtraEdges1[0]]+NodesCoordY[TempExtraEdges2[0]])/2;
        NrOfNewNodes = 1;
        
        for (i=1; i<(3*NrOfFirstZoneCandidateElements); i++){
            XCoordOfTempEdgeNode1 = NodesCoordX[TempExtraEdges1[i]];
            YCoordOfTempEdgeNode1 = NodesCoordY[TempExtraEdges1[i]];
            XCoordOfTempEdgeNode2 = NodesCoordX[TempExtraEdges2[i]];
            YCoordOfTempEdgeNode2 = NodesCoordY[TempExtraEdges2[i]];
            XTempNodeCoord = (XCoordOfTempEdgeNode1 + XCoordOfTempEdgeNode2)/2;
            YTempNodeCoord = (YCoordOfTempEdgeNode1 + YCoordOfTempEdgeNode2)/2;
            TempEdgeLength = SQRT((XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)*(XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)
            +(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1)*(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1));
            if (TempEdgeLength > Tol){
                // Check if new point is very close to already added new point
                AddCoordFlag = 1;
                for (j=0; j<NrOfNewNodes; j++){
                    DistNodes = SQRT((NewNodesCoordX[j]-XTempNodeCoord)*(NewNodesCoordX[j]-XTempNodeCoord)
                    + (NewNodesCoordY[j]-YTempNodeCoord)*(NewNodesCoordY[j]-YTempNodeCoord));
                    if (DistNodes < 2*EPSILON){
                        AddCoordFlag = 0;
                        break;
                    }
                }
                // If not then add new point
                if (AddCoordFlag == 1){
                    NewNodesCoordX[NrOfNewNodes] = XTempNodeCoord;
                    NewNodesCoordY[NrOfNewNodes] = YTempNodeCoord;
                    NrOfNewNodes = NrOfNewNodes+1;
                }
            }
        }
        
        // removing the first new node if the edge is too short
        XCoordOfTempEdgeNode1 = NodesCoordX[TempExtraEdges1[0]];
        YCoordOfTempEdgeNode1 = NodesCoordY[TempExtraEdges1[0]];
        XCoordOfTempEdgeNode2 = NodesCoordX[TempExtraEdges2[0]];
        YCoordOfTempEdgeNode2 = NodesCoordY[TempExtraEdges2[0]];
        TempEdgeLength = SQRT((XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)*(XCoordOfTempEdgeNode2-XCoordOfTempEdgeNode1)
        +(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1)*(YCoordOfTempEdgeNode2-YCoordOfTempEdgeNode1));
        if (TempEdgeLength<Tol){
            NrOfNewNodes = NrOfNewNodes-1;
            for (i=0; i<NrOfNewNodes; i++){
                NewNodesCoordX[i] = NewNodesCoordX[i+1];
                NewNodesCoordY[i] = NewNodesCoordY[i+1];
            }
        }
//         printf("NrOfNewNodes=%ld\n", NrOfNewNodes);
//     return;
        
        //  For triangles that showed up exactly once check if they are skinny        
        REAL TempLengths[3] = {0};
        REAL XNode1Coord, XNode2Coord, XNode3Coord;
        REAL YNode1Coord, YNode2Coord, YNode3Coord;
        REAL SkinnyTriangle = 3.0;
        for (i=0; i<NrOfElements; i++){
            if (FlagOfSecondZoneCandidateElements[i] == 1){
                XNode1Coord = NodesCoordX[Elements1[i]];
                YNode1Coord = NodesCoordY[Elements1[i]];
                XNode2Coord = NodesCoordX[Elements2[i]];
                YNode2Coord = NodesCoordY[Elements2[i]];
                XNode3Coord = NodesCoordX[Elements3[i]];
                YNode3Coord = NodesCoordY[Elements3[i]];
                
                TempLengths[0] = SQRT((XNode2Coord-XNode1Coord)*(XNode2Coord-XNode1Coord)
                + (YNode2Coord- YNode1Coord)*(YNode2Coord- YNode1Coord));
                TempLengths[1] = SQRT((XNode3Coord-XNode2Coord)*(XNode3Coord-XNode2Coord)
                + (YNode3Coord- YNode2Coord)*(YNode3Coord- YNode2Coord));
                TempLengths[2] = SQRT((XNode1Coord-XNode3Coord)*(XNode1Coord-XNode3Coord)
                + (YNode1Coord- YNode3Coord)*(YNode1Coord- YNode3Coord));
                // If skinny then add midpoint of triangle as newnode
                if (misc_max(3,&TempLengths[0])/misc_min(3,&TempLengths[0]) > SkinnyTriangle){
                    NewNodesCoordX[NrOfNewNodes] = (XNode1Coord+XNode2Coord+XNode3Coord)/3;
                    NewNodesCoordY[NrOfNewNodes] = (YNode1Coord+YNode2Coord+YNode3Coord)/3;
                    NrOfNewNodes = NrOfNewNodes+1;
                }
            }
        }
//         printf("NrOfNewNodes=%ld\n", NrOfNewNodes);

    }
//     printf("NrOfCandidateEdges=%ld\n",NrOfCandidateEdges );

    if (!converged_flag)
        WARN("Desired tolerance not achieved. Try increasing the number of iterations.");
    
// Evaluation of contour edges from all candidates edges
    NrOfCandidateElements = 0;
    ArrayOfCandidateElements= realloc(ArrayOfCandidateElements, 3*NrOfElements*sizeof(UINT));
    if (ArrayOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    ret_code = delaunay_triangulation_edgeAttachments(NrOfElements,Elements1, Elements2, Elements3,
            NrOfCandidateEdges, CandidateEdges1, CandidateEdges2, &NrOfCandidateElements, ArrayOfCandidateElements);
    if (ret_code != SUCCESS){
        ret_code = E_SUBROUTINE(ret_code);
        CHECK_RETCODE(ret_code, leave_fun);
    }
//    printf("NrOfCandidateElements=%ld\n",NrOfCandidateElements);
// Removing repeated mention of triangles
    ElementsCount = realloc(ElementsCount, NrOfElements*sizeof(UINT));
    if (ElementsCount == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i=0; i<NrOfElements; i++)
        ElementsCount[i] = 0;
    
    for (i=0; i<NrOfCandidateElements; i++)
        ElementsCount[ArrayOfCandidateElements[i]] = 1;
    
    IDOfCandidateElements = malloc(NrOfCandidateElements*sizeof(UINT));
    if (IDOfCandidateElements == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    NrOfCandidateElements = 0; // without reptition
    for (i=0; i<NrOfElements; i++){
        if (ElementsCount[i] == 1){
            IDOfCandidateElements[NrOfCandidateElements] = i;
            NrOfCandidateElements = NrOfCandidateElements+1;
            
        }
    }
//     printf("NrOfCandidateElements=%ld\n",NrOfCandidateElements);
//
// % Collecting all the edges of all candidate triangles
// TempEdges1=zeros(NoOfCandidateElements*3,1);
// TempEdges2=zeros(NoOfCandidateElements*3,1);
    TempExtraEdges1 = realloc(TempExtraEdges1, NrOfCandidateElements*3*sizeof(UINT));
    TempExtraEdges2 = realloc(TempExtraEdges2, NrOfCandidateElements*3*sizeof(UINT));
    if (TempExtraEdges1 == NULL || TempExtraEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    for (i=0; i<NrOfCandidateElements; i++){
        TempExtraEdges1[i*3] = Elements1[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3] = Elements2[IDOfCandidateElements[i]];
        TempExtraEdges1[i*3 + 1] = Elements2[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3 + 1] = Elements3[IDOfCandidateElements[i]];
        TempExtraEdges1[i*3 + 2] = Elements3[IDOfCandidateElements[i]];
        TempExtraEdges2[i*3 + 2] = Elements1[IDOfCandidateElements[i]];
    }
    
    // Removing edges to keep only outer contours
    
    UINT NrOfContourEdges = 0;
    MultiplicationOfTempEdges = calloc(3*NrOfCandidateElements, sizeof(UINT)); //NOTE calloc
    ContourEdges1 = malloc(3*NrOfCandidateElements*sizeof(UINT));
    ContourEdges2 = malloc(3*NrOfCandidateElements*sizeof(UINT));
    if (MultiplicationOfTempEdges == NULL || ContourEdges1 == NULL ||
            ContourEdges2 == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    INT NrOfEdge;
    // Possibly simplify this
    for (i=0; i<(3*NrOfCandidateElements); i++){
        if (MultiplicationOfTempEdges[i] == 0){
            NrOfEdge = -1;
            for (j=0; j<(3*NrOfCandidateElements); j++){
                if (j != i){
                    if ((TempExtraEdges2[j] == TempExtraEdges1[i] && TempExtraEdges1[j] == TempExtraEdges2[i])
                            || (TempExtraEdges1[j] == TempExtraEdges1[i] && TempExtraEdges2[j] == TempExtraEdges2[i])){
                        NrOfEdge = j;
                        MultiplicationOfTempEdges[NrOfEdge] = 2;
                    }
                }
            }
            if (NrOfEdge == -1){
                MultiplicationOfTempEdges[i] = 1;
                ContourEdges1[NrOfContourEdges] = TempExtraEdges1[i];
                ContourEdges2[NrOfContourEdges] = TempExtraEdges2[i];
                NrOfContourEdges = NrOfContourEdges+1;
            }else
                MultiplicationOfTempEdges[i] = 2;
        }
    }
    
//     printf("NrOfContourEdges=%ld\n",NrOfContourEdges);
    
    // evaluation of the regions
    NrOfNodesOfRegions = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    Regions = calloc(5*NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    IndexOfNextEdge = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    TempNodes = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    SideFlag = calloc(NrOfContourEdges, sizeof(UINT)); //NOTE calloc
    if (NrOfNodesOfRegions == NULL || Regions == NULL || IndexOfNextEdge == NULL
            || TempNodes == NULL || SideFlag == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    
    UINT NodesCount = 1;
    UINT NrOfRegions = 1;
    Regions[NodesCount-1] = ContourEdges1[0];
    NrOfNodesOfRegions[NrOfRegions-1]++;
    UINT RefNode = ContourEdges2[0];
    UINT ContourEdgesStart = 1;
    UINT NrOfIndexOfNextEdge;
    UINT PrevNode, Index;
    while (ContourEdgesStart < NrOfContourEdges){
        NrOfIndexOfNextEdge = 0;
        for (j=ContourEdgesStart; j<NrOfContourEdges; j++){
            if (ContourEdges1[j] == RefNode){
                IndexOfNextEdge[NrOfIndexOfNextEdge] = j;
                SideFlag[NrOfIndexOfNextEdge] = 1;
                NrOfIndexOfNextEdge = NrOfIndexOfNextEdge+1;
            }else if (ContourEdges2[j] == RefNode){
                IndexOfNextEdge[NrOfIndexOfNextEdge] = j;
                SideFlag[NrOfIndexOfNextEdge] = 2;
                NrOfIndexOfNextEdge = NrOfIndexOfNextEdge+1;
            }
        }
        if (NrOfIndexOfNextEdge == 0){
            NodesCount++;
            Regions[NodesCount-1] = RefNode;
            NrOfNodesOfRegions[NrOfRegions-1]++;
            
            if (ContourEdgesStart< NrOfContourEdges){
                NrOfRegions++;
                NodesCount++;
                Regions[NodesCount-1] = ContourEdges1[ContourEdgesStart];
                NrOfNodesOfRegions[NrOfRegions-1]++;
                RefNode = ContourEdges2[ContourEdgesStart];
                ContourEdgesStart = ContourEdgesStart+1;
            }
        }else{
            if (NrOfIndexOfNextEdge > 1){
                PrevNode = Regions[NodesCount-1];
                for (i=0; i<NrOfIndexOfNextEdge; i++){
                    if (SideFlag[i] == 1)
                        TempNodes[i] = ContourEdges2[IndexOfNextEdge[NrOfIndexOfNextEdge]];
                    else if (SideFlag[i] == 2)
                        TempNodes[i] = ContourEdges1[IndexOfNextEdge[NrOfIndexOfNextEdge]];
                }
                Index = FindNextNode(NodesCoordX, NodesCoordY, PrevNode, RefNode, TempNodes, NrOfIndexOfNextEdge);
                IndexOfNextEdge[0] = IndexOfNextEdge[Index];
                SideFlag[0] = SideFlag[Index];
                NrOfIndexOfNextEdge = 1;
            }
            NodesCount++;
            Regions[NodesCount-1] = RefNode;
            NrOfNodesOfRegions[NrOfRegions-1]++;
            if (SideFlag[0] == 1)
                RefNode = ContourEdges2[IndexOfNextEdge[0]];
            else if (SideFlag[0] == 2)
                RefNode = ContourEdges1[IndexOfNextEdge[0]];
            
            NrOfContourEdges = NrOfContourEdges - 1;
            for (j=ContourEdgesStart; j<NrOfContourEdges; j++){
                if (j>=IndexOfNextEdge[0]){
                    ContourEdges1[j] = ContourEdges1[j+1];
                    ContourEdges2[j] = ContourEdges2[j+1];
                }
            }
        }
    }
    NodesCount++;
    Regions[NodesCount-1] = RefNode;
    NrOfNodesOfRegions[NrOfRegions-1]++;
    //printf("NrOfRegions=%ld\n",NrOfRegions);
    //for (i=0; i<NrOfRegions; i++)
      // printf("Nodes in region %ld = %ld\n",i, NrOfNodesOfRegions[i]);
    
    z_m = calloc(NrOfRegions, sizeof(INT)); //NOTE calloc
    z = calloc(NrOfRegions, sizeof(COMPLEX)); //NOTE calloc
    contour_orientation = malloc(NrOfRegions*sizeof(INT));
    if (z_m == NULL || z == NULL || contour_orientation == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    NodesCount = 0;
    INT dQ;
    for (i=0; i<NrOfRegions; i++){
        for (j=1; j<NrOfNodesOfRegions[i]; j++){
            dQ = Quadrants[Regions[NodesCount+j]] - Quadrants[Regions[NodesCount+j-1]];
            //printf("dQ=%d\n",z_m[i]);
            if (dQ == 3)
                z_m[i] = z_m[i] - 1;
            else if (dQ == -3)
                z_m[i] = z_m[i] + 1;
            else if (ABS(dQ) == 2){
                z_m[i] = 0;
                break;
            }
            else
                z_m[i] = z_m[i] + dQ;
            
        }

        z_m[i] = z_m[i]/4;
        //printf("z_m[i]=%d\n",z_m[i]);
        z[i] = 0;
        contour_orientation[i] = 0;
        // If z_m[i]==0 then discard region
        if (z_m[i] != 0){
           for (j=0; j<NrOfNodesOfRegions[i]; j++){
               z[i] = z[i]+(NodesCoordX[Regions[NodesCount+j]] + I*NodesCoordY[Regions[NodesCount+j]]);
               //printf("x=%1.5e,y=%1.5e\n",NodesCoordX[Regions[NodesCount+j]],NodesCoordY[Regions[NodesCount+j]]);
           }

           z[i] = z[i]/NrOfNodesOfRegions[i]; // Approximation of zero or pole

           // We now figure out the orientation of the contour
           // We do this by setting computed z[i] as the origin and computing the qudrants of the points
           // in the current region. By looking at the order of the quadrants we get the contour orientation.

           free(ComplexPoints);
           ComplexPoints = malloc(NrOfNodesOfRegions[i] * sizeof(COMPLEX)); //TODO : Reuse previously allocated memory
           free(NewQuadrants); 
           NewQuadrants = malloc(NrOfNodesOfRegions[i] * sizeof(UINT));
           free(NewAngles);
           NewAngles = malloc(NrOfNodesOfRegions[i] * sizeof(REAL));
           if (NewQuadrants == NULL || ComplexPoints == NULL || NewAngles == NULL){
              ret_code = E_NOMEM;
              goto leave_fun;
           }

           for (j=0; j<NrOfNodesOfRegions[i]; j++)
               ComplexPoints[j] = (NodesCoordX[Regions[NodesCount+j]] + I*NodesCoordY[Regions[NodesCount+j]]) - z[i];

           ret_code = misc_quadrant(NrOfNodesOfRegions[i], ComplexPoints,
                NewQuadrants, NewAngles);
           if (ret_code != SUCCESS){
              ret_code = E_SUBROUTINE(ret_code);
              CHECK_RETCODE(ret_code, leave_fun);
           }

           // Clockwise contour with negative z_m implies a zero
           // Anti-clockwise contour with negative z_m implies a pole
           // Clockwise contour with positive z_m implies a pole
           // Anti-clockwise contour with positive z_m implies a zero


           INT Qj, Qjp1;
           for (j=0; j<NrOfNodesOfRegions[i]-1; j++){
               Qjp1 = NewQuadrants[j+1];
               Qj = NewQuadrants[j];
               if (ABS(Qjp1 - Qj) == 1){
                  if ((NewQuadrants[j+1] == 1 && NewQuadrants[j] == 4) || NewQuadrants[j+1]>NewQuadrants[j]){
                     contour_orientation[i] = -1; // Anti-clockwise contour
                     break; 
                  }
                  else{
                     contour_orientation[i] = 1;  // Clockwise contour
                     break; 
                  }
               }
           } 


        }
        NodesCount = NodesCount+NrOfNodesOfRegions[i];
    }
    
    if (roots != NULL){
    	j = 0;
    	for (i=0; i<NrOfRegions; i++){
    	    if (*K_roots_ptr>j){
      	        if (z_m[i]*contour_orientation[i] < 0){ // Zeros
         	    roots[j] = z[i];  
                    if (roots_multiplicity != NULL)
      	                roots_multiplicity[j] = ABS(z_m[i]);
                    j++;
                }
            }  
            else{
                WARN("Ran out of memory for roots.");
                break;
            } 
        }
        *K_roots_ptr = j;
    }


    if (poles != NULL){    
        j = 0;
        for (i=0; i<NrOfRegions; i++){
            if (*K_poles_ptr>j){
               if (z_m[i]*contour_orientation[i] > 0){ // Poles
                  poles[j] = z[i];
                  if (poles_multiplicity != NULL)
                     poles_multiplicity[j] = ABS(z_m[i]);
                  j++;
               }
            }
            else{
                WARN("Ran out of memory for poles.");
                break;
            }
        }
        *K_poles_ptr = j;
    }

   
    leave_fun:
        free(NewNodesCoordX);
        free(NewNodesCoordY);
        free(NodesCoordX);
        free(NodesCoordY);
        free(FuntionValues);
        free(ComplexNodes);
        free(Quadrants);
        free(Angles);
        free(ValidEdges);
        free(ValidEdges1);
        free(ValidEdges2);
        free(CandidateEdges1);
        free(CandidateEdges2);
        free(Elements1);
        free(Elements2);
        free(Elements3);
        free(FlagOfFirstZoneCandidateElements);
        free(FlagOfSecondZoneCandidateElements);
        free(ElementsCount);
        free(IDOfCandidateElements);
        free(TempExtraEdges1);
        free(TempExtraEdges2);
        free(ContourEdges1);
        free(ContourEdges2);
        free(MultiplicationOfTempEdges);
        free(SideFlag);
        free(TempNodes);
        free(IndexOfNextEdge);
        free(Regions);
        free(NrOfNodesOfRegions);
        free(z);
        free(z_m);
        free(ArrayOfCandidateElements);
        free(CandidateNodes);
        free(Temp);
        free(ComplexPoints);
        free(NewQuadrants);
        free(NewAngles);
        free(contour_orientation);
        if (DT_built == 1){
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
        }
        return ret_code;

}


