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

#include "fnft__delaunay_triangulation.h"

INT CircumCircle(REAL *X, REAL *Y, REAL *xcoord, REAL *ycoord, REAL *R){
  
    REAL s1 = SQRT((X[1]-X[0])*(X[1]-X[0]) + (Y[1]-Y[0])*(Y[1]-Y[0]));
    REAL s2 = SQRT((X[2]-X[1])*(X[2]-X[1]) + (Y[2]-Y[1])*(Y[2]-Y[1]));
    REAL s3 = SQRT((X[0]-X[2])*(X[0]-X[2]) + (Y[0]-Y[2])*(Y[0]-Y[2]));
    
    REAL denom = (s1+s2+s3)*(-s1+s2+s3)*(s1-s2+s3)*(s1+s2-s3);
    if (denom < 0.0) 
        return E_ASSERTION_FAILED;
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    R[0] = (s1*s2*s3)/SQRT(denom);

    denom = 2*s1*s3;
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    REAL ang1 = ACOS((s1*s1+s3*s3-s2*s2)/denom);

    denom = 2*s1*s2;
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    REAL ang2 = ACOS((s1*s1+s2*s2-s3*s3)/denom);

    denom = 2*s2*s3;
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    REAL ang3 = ACOS((s2*s2+s3*s3-s1*s1)/denom);
   
    denom = SIN(2*ang1)+SIN(2*ang2)+SIN(2*ang3);
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    xcoord[0] = (X[0]*SIN(2*ang1) + X[1]*SIN(2*ang2) + X[2]*SIN(2*ang3))/denom;

    denom = SIN(2*ang1)+SIN(2*ang2)+SIN(2*ang3);
    if (denom == 0.0)
        return E_DIV_BY_ZERO;
    ycoord[0] = (Y[0]*SIN(2*ang1) + Y[1]*SIN(2*ang2) + Y[2]*SIN(2*ang3))/denom;

    return SUCCESS;
}

INT delaunay_triangulation(fnft__delaunay_triangulation_data_t * DT_ptr,
        UINT NodesMax, UINT NrOfNewNodes,
        REAL * NewNodesCoordX, REAL * NewNodesCoordY, UINT * status_flag_ptr)
{
    
    if (DT_ptr->NrOfNodes != 0){
        if (DT_ptr-> NodesMax != NodesMax)
            return E_INVALID_ARGUMENT(NodesMax);
    }
    if (*status_flag_ptr > 0)
        return E_INVALID_ARGUMENT(status_flag_ptr);
    if (NrOfNewNodes == 0)
        return SUCCESS;
    if (NewNodesCoordX == NULL)
        return E_INVALID_ARGUMENT(NewNodesCoordX);
    if (NewNodesCoordY == NULL)
        return E_INVALID_ARGUMENT(NewNodesCoordY);
    
    
    INT ret_code = SUCCESS;
    UINT i, j, add_edge_flag, NewEdge1,NewEdge2;
    UINT *IdxOfBadTriangles = NULL, *PolygonEdges1 = NULL,
            *PolygonEdges2 = NULL, *RepeatedEdge = NULL;
    REAL Xtmp[3], Ytmp[3], CCX, CCY, CCR, m;
    
    if (DT_ptr->NodesMax == 0){
        
        // Allocating memory and initializing certain values
        DT_ptr->NodesMax = NodesMax;
        REAL ye = misc_max(NrOfNewNodes, NewNodesCoordY);
        REAL yb = misc_min(NrOfNewNodes, NewNodesCoordY);
        REAL xb = misc_min(NrOfNewNodes, NewNodesCoordX);
        REAL xe = misc_max(NrOfNewNodes, NewNodesCoordX);
        REAL xmax = xe-xb;
        REAL ymax = ye-yb;     
                
        DT_ptr->NodeOfTriangles1 = malloc(2*NodesMax * sizeof(UINT));
        DT_ptr->NodeOfTriangles2 = malloc(2*NodesMax * sizeof(UINT));
        DT_ptr->NodeOfTriangles3 = malloc(2*NodesMax * sizeof(UINT));
        if (DT_ptr->NodeOfTriangles1 == NULL || DT_ptr->NodeOfTriangles2 == NULL
                || DT_ptr->NodeOfTriangles3 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        DT_ptr->EdgesOfTriangles1 = malloc(2*NodesMax * sizeof(UINT));
        DT_ptr->EdgesOfTriangles2 = malloc(2*NodesMax * sizeof(UINT));
        DT_ptr->EdgesOfTriangles3 = malloc(2*NodesMax * sizeof(UINT));
        if (DT_ptr->EdgesOfTriangles1 == NULL || DT_ptr->EdgesOfTriangles2 == NULL
                || DT_ptr->EdgesOfTriangles3 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        DT_ptr->Edges1 = malloc(4*NodesMax * sizeof(UINT));
        DT_ptr->Edges2 = malloc(4*NodesMax * sizeof(UINT));
        if (DT_ptr->Edges1 == NULL || DT_ptr->Edges2 == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        DT_ptr->NodesCoordX = malloc(NodesMax * sizeof(REAL));
        DT_ptr->NodesCoordY = malloc(NodesMax * sizeof(REAL));
        if (DT_ptr->NodesCoordX == NULL || DT_ptr->NodesCoordY == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }

        //    Super rectangle
        DT_ptr->NodesCoordX[0] = xb-xmax/2;
        DT_ptr->NodesCoordY[0] = yb-ymax/2;
        DT_ptr->NodesCoordX[1] = xe+xmax/2;
        DT_ptr->NodesCoordY[1] = yb-ymax/2;
        DT_ptr->NodesCoordX[2] = xb-xmax/2;
        DT_ptr->NodesCoordY[2] = ye+ymax/2;
        DT_ptr->NodesCoordX[3] = xe+xmax/2;
        DT_ptr->NodesCoordY[3] = ye+ymax/2;
        DT_ptr->NrOfNodes = 4;
        
        DT_ptr->NrOfTriangles = 2;
        DT_ptr->NodeOfTriangles1[0] = 0;//First node of first triangle
        DT_ptr->NodeOfTriangles2[0] = 1;//Second node of first triangle
        DT_ptr->NodeOfTriangles3[0] = 2;
        DT_ptr->NodeOfTriangles1[1] = 1;//First node of second triangle
        DT_ptr->NodeOfTriangles2[1] = 2;//Second node of second triangle
        DT_ptr->NodeOfTriangles3[1] = 3;
        DT_ptr->NrOfNodesToRemove = (4);
        //     % NOTE: Edges1(k)<Edges2(k)
        DT_ptr->Edges1[0] = 0; //Nodes
        DT_ptr->Edges2[0] = 1; //Nodes
        DT_ptr->Edges1[1] = 1; //Nodes
        DT_ptr->Edges2[1] = 2; //Nodes
        DT_ptr->Edges1[2] = 0; //Nodes
        DT_ptr->Edges2[2] = 2; //Nodes
        DT_ptr->Edges1[3] = 2; //Nodes
        DT_ptr->Edges2[3] = 3; //Nodes
        DT_ptr->Edges1[4] = 1; //Nodes
        DT_ptr->Edges2[4] = 3; //Nodes
        DT_ptr->NrOfEdges = 5;
        
        DT_ptr->EdgesOfTriangles1[0] = 0; //Edge numbers
        DT_ptr->EdgesOfTriangles2[0] = 1;
        DT_ptr->EdgesOfTriangles3[0] = 2;
        DT_ptr->EdgesOfTriangles1[1] = 1; //Edge numbers
        DT_ptr->EdgesOfTriangles2[1] = 3;
        DT_ptr->EdgesOfTriangles3[1] = 4;
        
        DT_ptr->XCoordOfCCOfTriangles = malloc(2*NodesMax * sizeof(REAL));
        DT_ptr->YCoordOfCCOfTriangles = malloc(2*NodesMax * sizeof(REAL));
        DT_ptr->RadiusOfCCOfTriangles = malloc(2*NodesMax * sizeof(REAL));
        if (DT_ptr->XCoordOfCCOfTriangles == NULL || DT_ptr->YCoordOfCCOfTriangles == NULL
                || DT_ptr->RadiusOfCCOfTriangles == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        
        for (i=0; i<DT_ptr->NrOfTriangles; i++){
            Xtmp[0] = DT_ptr->NodesCoordX[DT_ptr->NodeOfTriangles1[i]];
            Xtmp[1] = DT_ptr->NodesCoordX[DT_ptr->NodeOfTriangles2[i]];
            Xtmp[2] = DT_ptr->NodesCoordX[DT_ptr->NodeOfTriangles3[i]];
            Ytmp[0] = DT_ptr->NodesCoordY[DT_ptr->NodeOfTriangles1[i]];
            Ytmp[1] = DT_ptr->NodesCoordY[DT_ptr->NodeOfTriangles2[i]];
            Ytmp[2] = DT_ptr->NodesCoordY[DT_ptr->NodeOfTriangles3[i]];
            
            ret_code = CircumCircle(&Xtmp[0], &Ytmp[0], &CCX, &CCY, &CCR);
            CHECK_RETCODE(ret_code, leave_fun); 
            DT_ptr->XCoordOfCCOfTriangles[i] = CCX;
            DT_ptr->YCoordOfCCOfTriangles[i] = CCY;
            DT_ptr->RadiusOfCCOfTriangles[i] = CCR;
        }
        
        DT_ptr->StatusOfTriangles = malloc(2*NodesMax * sizeof(UINT));
        DT_ptr->CheckEdge = malloc(2*NodesMax * sizeof(UINT));
        if (DT_ptr->StatusOfTriangles == NULL || DT_ptr->CheckEdge == NULL){
            ret_code = E_NOMEM;
            goto leave_fun;
        }
        DT_ptr->StatusOfTriangles[0] = 1; // 1 means good
        DT_ptr->StatusOfTriangles[1] = 1; // 1 means good
        
        
        DT_ptr->CheckEdge[0] = 0; // 1 means yet to add/assign edge to triangle
        DT_ptr->CheckEdge[1] = 0;
        
        DT_ptr->RemoveSuperRectangle = 1;
        
    }
    
    if (DT_ptr->NodesMax != NodesMax){
        ret_code = E_ASSERTION_FAILED;
        goto leave_fun;
    }

    REAL * NodesCoordX = DT_ptr->NodesCoordX;
    REAL * NodesCoordY = DT_ptr->NodesCoordY;
    UINT * NodeOfTriangles1 = DT_ptr->NodeOfTriangles1;
    UINT * NodeOfTriangles2 = DT_ptr->NodeOfTriangles2;
    UINT * NodeOfTriangles3 = DT_ptr->NodeOfTriangles3;
    UINT * Edges1 = DT_ptr->Edges1;
    UINT * Edges2 = DT_ptr->Edges2;
    UINT * CheckEdge = DT_ptr->CheckEdge;
    
    REAL * XCoordOfCCOfTriangles = DT_ptr->XCoordOfCCOfTriangles;
    REAL * YCoordOfCCOfTriangles = DT_ptr->YCoordOfCCOfTriangles;
    REAL * RadiusOfCCOfTriangles = DT_ptr->RadiusOfCCOfTriangles;
    
    UINT * StatusOfTriangles = DT_ptr->StatusOfTriangles;
    UINT NrOfEdges = DT_ptr->NrOfEdges;
    UINT * EdgesOfTriangles1 = DT_ptr->EdgesOfTriangles1;
    UINT * EdgesOfTriangles2 = DT_ptr->EdgesOfTriangles2;
    UINT * EdgesOfTriangles3 = DT_ptr->EdgesOfTriangles3;
    UINT NrOfNodes = DT_ptr->NrOfNodes;
    UINT NrOfTriangles = DT_ptr->NrOfTriangles;
    UINT NrOfNodesToRemove = DT_ptr->NrOfNodesToRemove;
    
    IdxOfBadTriangles = malloc(2*NodesMax * sizeof(UINT));
    PolygonEdges1 = malloc(2*NodesMax * sizeof(UINT));
    PolygonEdges2 = malloc(2*NodesMax * sizeof(UINT));
    RepeatedEdge = malloc(2*NodesMax * sizeof(UINT));
    if (IdxOfBadTriangles == NULL || PolygonEdges1 == NULL
            || PolygonEdges2 == NULL || RepeatedEdge == NULL){
        ret_code = E_NOMEM;
        goto leave_fun;
    }
    
    UINT k, new_node, NrOfBadTriangles, IdxOfTriangle, NrOfPolygonEdges, idxi;
    REAL dist;
    for (k=0; k<NrOfNewNodes; k++){
        
        new_node = 1;
        // Check if new node is very close to old node
        for (i=0; i<NrOfNodes; i++){
            if ((FABS(NodesCoordX[i]-NewNodesCoordX[k])<2*EPSILON) && (FABS(NodesCoordY[i]-NewNodesCoordY[k])<2*EPSILON)){
                new_node = 0;
//             warning("Removing duplicate point");
                *status_flag_ptr = 1;
                ret_code = SUCCESS;
                goto leave_fun;
            }
        }
        
        if (NrOfTriangles > (2*NodesMax)-10){
            *status_flag_ptr = 2;
            ret_code = SUCCESS;
            goto leave_fun;
        }
        if (NrOfEdges > (4*NodesMax)-10){
            *status_flag_ptr = 3;
            ret_code = SUCCESS;
            goto leave_fun;
        }
        if (new_node == 1){
            // Find triangles inside whose circumcircle the new point lies
            NrOfBadTriangles = 0;
            for  (IdxOfTriangle=0; IdxOfTriangle<NrOfTriangles; IdxOfTriangle++){
                if (StatusOfTriangles[IdxOfTriangle] == 1){
                    dist = ((XCoordOfCCOfTriangles[IdxOfTriangle]-NewNodesCoordX[k])*(XCoordOfCCOfTriangles[IdxOfTriangle]-NewNodesCoordX[k])
                    +(YCoordOfCCOfTriangles[IdxOfTriangle]-NewNodesCoordY[k])*(YCoordOfCCOfTriangles[IdxOfTriangle]-NewNodesCoordY[k]));
                    
                    if (dist < RadiusOfCCOfTriangles[IdxOfTriangle]*RadiusOfCCOfTriangles[IdxOfTriangle]){
                        IdxOfBadTriangles[NrOfBadTriangles] = IdxOfTriangle;
                        NrOfBadTriangles = NrOfBadTriangles+1;
                    }
                }
            }
            if (NrOfBadTriangles == 0){
                ret_code = E_ASSERTION_FAILED;
                goto leave_fun;
            }
            // List all edges of bad triangles
            NrOfPolygonEdges = 0;
            for (i=0; i<NrOfBadTriangles; i++){
                idxi = IdxOfBadTriangles[i];
                // First edge
                if (NodeOfTriangles1[idxi] < NodeOfTriangles2[idxi]){
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles1[idxi];
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles2[idxi];
                }else{
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles1[idxi];
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles2[idxi];
                }
                NrOfPolygonEdges = NrOfPolygonEdges+1;
                
                // Second edge
                if (NodeOfTriangles2[idxi] < NodeOfTriangles3[idxi]){
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles2[idxi];
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles3[idxi];
                }else{
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles2[idxi];
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles3[idxi];
                }
                NrOfPolygonEdges = NrOfPolygonEdges+1;
                
                // Third edge
                if (NodeOfTriangles3[idxi] < NodeOfTriangles1[idxi]){
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles3[idxi];
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles1[idxi];
                }else{
                    PolygonEdges2[NrOfPolygonEdges] = NodeOfTriangles3[idxi];
                    PolygonEdges1[NrOfPolygonEdges] = NodeOfTriangles1[idxi];
                }
                NrOfPolygonEdges = NrOfPolygonEdges+1;
            }
            
            // Setting all repeat flags to 0
            for (i=0; i<NrOfPolygonEdges; i++)
                RepeatedEdge[i] = 0;
            
            // Mark double edges
            for (i=0; i<NrOfPolygonEdges; i++){
                for (j=0; j<NrOfPolygonEdges; j++){
                    if (i != j){
                        if (PolygonEdges1[i] == PolygonEdges1[j]){
                            if (PolygonEdges2[i] == PolygonEdges2[j]){
                                RepeatedEdge[i] = 1;
                                break;
                            }
                        }
                    }
                }
            }
            
            j = 0;
            for (i=0; i<NrOfPolygonEdges; i++){
                if (RepeatedEdge[i] == 0){
                    PolygonEdges1[j] = PolygonEdges1[i];
                    PolygonEdges2[j] = PolygonEdges2[i];
                    j = j+1;
                }
            }
            NrOfPolygonEdges = j;
            
            // Set status of bad triangles to 0
            for (i=0; i<NrOfBadTriangles; i++)
                StatusOfTriangles[IdxOfBadTriangles[i]] = 0;
            
            
            if (NrOfPolygonEdges != 0){
                // Add new point as node
                NodesCoordX[NrOfNodes] = NewNodesCoordX[k];
                NodesCoordY[NrOfNodes] = NewNodesCoordY[k];
                NrOfNodes = NrOfNodes+1;
            }
            
            // Checking if new node lies on any edge
            // For everyother edge add new triangle
            for (i=0; i<NrOfPolygonEdges; i++){
                const REAL denom = NodesCoordX[PolygonEdges2[i]] - NodesCoordX[PolygonEdges1[i]]; 
                if (denom == 0.0) {
                    ret_code = E_DIV_BY_ZERO;
                    goto leave_fun;
                }
                m = (NodesCoordY[PolygonEdges2[i]] - NodesCoordY[PolygonEdges1[i]])/denom;
                if (FABS(NodesCoordY[NrOfNodes-1] - (m*(NodesCoordX[NrOfNodes-1] - NodesCoordX[PolygonEdges1[i]]) + NodesCoordY[PolygonEdges1[i]]))>2*EPSILON){
                    StatusOfTriangles[NrOfTriangles] = 1;
                    NodeOfTriangles1[NrOfTriangles] = PolygonEdges1[i];
                    NodeOfTriangles2[NrOfTriangles] = NrOfNodes-1; //new point is a vertex
                    NodeOfTriangles3[NrOfTriangles] = PolygonEdges2[i];
                    
                    // Calculate CC and R for new triangles
                    Xtmp[0] = NodesCoordX[NodeOfTriangles1[NrOfTriangles]];
                    Xtmp[1] = NodesCoordX[NodeOfTriangles2[NrOfTriangles]];
                    Xtmp[2] = NodesCoordX[NodeOfTriangles3[NrOfTriangles]];
                    Ytmp[0] = NodesCoordY[NodeOfTriangles1[NrOfTriangles]];
                    Ytmp[1] = NodesCoordY[NodeOfTriangles2[NrOfTriangles]];
                    Ytmp[2] = NodesCoordY[NodeOfTriangles3[NrOfTriangles]];
                    
                    ret_code = CircumCircle(&Xtmp[0], &Ytmp[0], &CCX, &CCY, &CCR);
                    CHECK_RETCODE(ret_code, leave_fun);
                    XCoordOfCCOfTriangles[NrOfTriangles] = CCX;
                    YCoordOfCCOfTriangles[NrOfTriangles] = CCY;
                    RadiusOfCCOfTriangles[NrOfTriangles] = CCR;
                    
                    if (CCR == INFINITY){
                        ret_code = E_ASSERTION_FAILED;
                        goto leave_fun;
                    }
                    CheckEdge[NrOfTriangles] = 1;
                    
                    NrOfTriangles = NrOfTriangles +1;
                }
            }
        }
    }
    //Remove super rectangle
    if (DT_ptr->RemoveSuperRectangle == 1){
        DT_ptr->RemoveSuperRectangle = 0;
        for (i=0; i<NrOfTriangles; i++){
            if (StatusOfTriangles[i] == 1){
                if (NodeOfTriangles1[i] < NrOfNodesToRemove
                        || NodeOfTriangles2[i] < NrOfNodesToRemove
                        || NodeOfTriangles3[i] < NrOfNodesToRemove)
                    StatusOfTriangles[i] = 0;
            }
        }
    }    
//         printf("NrOfEdges=%ld before adding new edges\n",NrOfEdges);  

    for (i=0; i<NrOfTriangles; i++){
        if (StatusOfTriangles[i] == 1 && CheckEdge[i] == 1){
            CheckEdge[i] = 0;
            // First edge
            if (NodeOfTriangles1[i] < NodeOfTriangles2[i]){
                NewEdge1 = NodeOfTriangles1[i];
                NewEdge2 = NodeOfTriangles2[i];
            }else{
                NewEdge2 = NodeOfTriangles1[i];
                NewEdge1 = NodeOfTriangles2[i];
            }
            add_edge_flag = 1;
            for (j=0; j<NrOfEdges; j++){
                if (NewEdge1 == Edges1[j]){
                    if (NewEdge2 == Edges2[j]){
                        add_edge_flag = 0;
                        EdgesOfTriangles1[i] = j;
                        break;
                    }
                }
            }
            if (add_edge_flag == 1){
                Edges1[NrOfEdges] = NewEdge1;
                Edges2[NrOfEdges] = NewEdge2;
                EdgesOfTriangles1[i] = NrOfEdges;
                NrOfEdges = NrOfEdges+1;
                
            }
            
            // Second edge
            if (NodeOfTriangles2[i] < NodeOfTriangles3[i]){
                NewEdge1 = NodeOfTriangles2[i];
                NewEdge2 = NodeOfTriangles3[i];
            }else{
                NewEdge2 = NodeOfTriangles2[i];
                NewEdge1 = NodeOfTriangles3[i];
            }
            add_edge_flag = 1;
            for (j=0; j<NrOfEdges; j++){
                if (NewEdge1 == Edges1[j]){
                    if (NewEdge2 == Edges2[j]){
                        add_edge_flag = 0;
                        EdgesOfTriangles2[i] = j;
                        break;
                    }
                }
            }
            if (add_edge_flag == 1){
                Edges1[NrOfEdges] = NewEdge1;
                Edges2[NrOfEdges] = NewEdge2;
                EdgesOfTriangles2[i] = NrOfEdges;
                NrOfEdges = NrOfEdges+1;
            }
            
            //Third edge
            if (NodeOfTriangles3[i] < NodeOfTriangles1[i]){
                NewEdge1 = NodeOfTriangles3[i];
                NewEdge2 = NodeOfTriangles1[i];
            }else{
                NewEdge2 = NodeOfTriangles3[i];
                NewEdge1 = NodeOfTriangles1[i];
            }
            add_edge_flag = 1;
            for (j=0; j<NrOfEdges; j++){
                if (NewEdge1 == Edges1[j]){
                    if (NewEdge2 == Edges2[j]){
                        add_edge_flag = 0;
                        EdgesOfTriangles3[i] = j;
                        break;
                    }
                }
            }
            if (add_edge_flag == 1){
                Edges1[NrOfEdges] = NewEdge1;
                Edges2[NrOfEdges] = NewEdge2;
                EdgesOfTriangles3[i] = NrOfEdges;
                NrOfEdges = NrOfEdges+1;
            }
            
        }
    }  

    // Cleaning up to save memory
    j = 0;
    for (i=0; i<NrOfTriangles; i++){
        if (StatusOfTriangles[i] == 1){
            NodeOfTriangles1[j] = NodeOfTriangles1[i];
            NodeOfTriangles2[j] = NodeOfTriangles2[i];
            NodeOfTriangles3[j] = NodeOfTriangles3[i];
            EdgesOfTriangles1[j] = EdgesOfTriangles1[i];
            EdgesOfTriangles2[j] = EdgesOfTriangles2[i];
            EdgesOfTriangles3[j] = EdgesOfTriangles3[i];
            XCoordOfCCOfTriangles[j] = XCoordOfCCOfTriangles[i];
            YCoordOfCCOfTriangles[j] = YCoordOfCCOfTriangles[i];
            RadiusOfCCOfTriangles[j] = RadiusOfCCOfTriangles[i];
            CheckEdge[j] = CheckEdge[i];
            j = j+1;
        }
    }
    NrOfTriangles = j;
    for (i=0; i<NrOfTriangles; i++)
        StatusOfTriangles[i] = 1;
    
    DT_ptr->NrOfEdges = NrOfEdges;
    DT_ptr->NrOfNodes = NrOfNodes;
    DT_ptr->NrOfTriangles = NrOfTriangles;
    DT_ptr->NrOfNodesToRemove = NrOfNodesToRemove;  
    
//     printf("NrOfNodes=%ld \n",NrOfNodes);
//     printf("NodesMax=%ld \n",NodesMax);
//     printf("NrOfTriangles=%ld \n",NrOfTriangles);
//     printf("NrOfEdges=%ld \n\n",NrOfEdges);  
    
    leave_fun:
        free(IdxOfBadTriangles);
        free(PolygonEdges1);
        free(PolygonEdges2);
        free(RepeatedEdge);
        return ret_code;
}

INT delaunay_triangulation_vertexAttachments(UINT NrOfElements,
        UINT *Elements1, UINT *Elements2, UINT *Elements3,
        UINT NrOfNodes, UINT *Nodes,
        UINT *NrOfCandidateElements_ptr,
        UINT *ArrayOfCandidateElements)
{
    if(NrOfElements == 0)
        return E_INVALID_ARGUMENT(NrOfElements);
    if (Elements1 == NULL)
        return E_INVALID_ARGUMENT(Elements1);
    if (Elements2 == NULL)
        return E_INVALID_ARGUMENT(Elements2);
    if (Elements3 == NULL)
        return E_INVALID_ARGUMENT(Elements3);
     if(NrOfNodes == 0)
        return E_INVALID_ARGUMENT(NrOfNodes);
    if (Nodes == NULL)
        return E_INVALID_ARGUMENT(Nodes);
    if (NrOfCandidateElements_ptr == NULL)
        return E_INVALID_ARGUMENT(NrOfCandidateElements_ptr);
    if (ArrayOfCandidateElements == NULL)
        return E_INVALID_ARGUMENT(ArrayOfCandidateElements);
    
    INT ret_code = SUCCESS;
    UINT count = 0, i, j;    
    for (i=0; i<NrOfNodes; i++){
        for (j=0; j<NrOfElements; j++){
            if (Elements1[j] == Nodes[i] || Elements2[j] == Nodes[i] || Elements3[j] == Nodes[i]){
                ArrayOfCandidateElements[count] = j;
                count = count+1;
            }
        }
    }
    *NrOfCandidateElements_ptr = count;    
    return ret_code;
}

INT delaunay_triangulation_edgeAttachments(UINT NrOfElements,
        UINT *Elements1, UINT *Elements2, UINT *Elements3,
        UINT NrOfNodes, UINT *Nodes1, UINT *Nodes2,
        UINT *NrOfCandidateElements_ptr,
        UINT *ArrayOfCandidateElements)
{
    if(NrOfElements == 0)
        return E_INVALID_ARGUMENT(NrOfElements);
    if (Elements1 == NULL)
        return E_INVALID_ARGUMENT(Elements1);
    if (Elements2 == NULL)
        return E_INVALID_ARGUMENT(Elements2);
    if (Elements3 == NULL)
        return E_INVALID_ARGUMENT(Elements3);
     if(NrOfNodes == 0)
        return E_INVALID_ARGUMENT(NrOfNodes);
    if (Nodes1 == NULL)
        return E_INVALID_ARGUMENT(Nodes1);
    if (Nodes2 == NULL)
        return E_INVALID_ARGUMENT(Nodes2);
    if (NrOfCandidateElements_ptr == NULL)
        return E_INVALID_ARGUMENT(NrOfCandidateElements_ptr);
    if (ArrayOfCandidateElements == NULL)
        return E_INVALID_ARGUMENT(ArrayOfCandidateElements);
    
    INT ret_code = SUCCESS;
    UINT count = 0, i, j;
    for (i=0; i<NrOfNodes; i++){
        for (j=0; j<NrOfElements; j++){
            if (Elements1[j] == Nodes1[i] || Elements2[j] == Nodes1[i] || Elements3[j] == Nodes1[i]){
                if (Elements1[j] == Nodes2[i] || Elements2[j] == Nodes2[i] || Elements3[j] == Nodes2[i]){
                    ArrayOfCandidateElements[count] = j;
                    count = count+1;
                }
            }
        }
    }
    *NrOfCandidateElements_ptr = count;    
    return ret_code;
}

INT delaunay_triangulation_rect_dom(REAL xb, REAL xe, REAL yb, REAL ye, REAL r,
        UINT * const NrOfNewNodesCoord_ptr, REAL ** NewNodesCoordX_ptr, REAL ** NewNodesCoordY_ptr)
{

    if (xe < xb)
        return E_INVALID_ARGUMENT(xe);
    if (ye < yb)
        return E_INVALID_ARGUMENT(ye);
    if (NewNodesCoordX_ptr == NULL)
        return E_INVALID_ARGUMENT(NewNodesCoordX_ptr);
    if (NewNodesCoordY_ptr == NULL)
        return E_INVALID_ARGUMENT(NewNodesCoordY_ptr);
    if (NrOfNewNodesCoord_ptr == NULL)
        return E_INVALID_ARGUMENT(NrOfNewNodesCoord_ptr);
    
    UINT i, j, NrOfNewNodesCoord;
    
    REAL X = xe-xb;
    REAL Y = ye-yb;
    
    UINT n = CEIL(Y/r+1);
    REAL dy = Y/(n-1);
    UINT m = CEIL(X/SQRT(r*r-dy*dy/4) +1);
    REAL dx = X/(m-1);
    
//     printf("m = %ld\n",m);
//     printf("n = %ld\n",n);
    UINT len = m*n+(UINT)FLOOR(m/2);
//     printf("len = %ld\n",len); 
    REAL *NewNodesCoordX = malloc(len * sizeof(REAL));
    REAL *NewNodesCoordY = malloc(len * sizeof(REAL));
    if (NewNodesCoordX == NULL || NewNodesCoordY == NULL)
        return E_NOMEM;
    
    NrOfNewNodesCoord = 0;
    for (i=0; i<m; i++){
        for (j=0; j<n; j++){
//             printf("NrOfNewNodesCoord +j = %ld    ",NrOfNewNodesCoord +j);
            NewNodesCoordX[NrOfNewNodesCoord +j] = xb + i*dx;
            NewNodesCoordY[NrOfNewNodesCoord +j] = yb + j*dy;
            if ( (j != n-1) && ((i+1)%2 == 0))
                NewNodesCoordY[NrOfNewNodesCoord+j] = NewNodesCoordY[NrOfNewNodesCoord+j] + 0.5*dy;
        }
        NrOfNewNodesCoord += n;
    }
//      printf("NrOfNewNodesCoord = %ld\n",NrOfNewNodesCoord); 
    
    for (i=1; i<m; i+=2){
        NewNodesCoordX[NrOfNewNodesCoord] = i*dx + xb;
        NewNodesCoordY[NrOfNewNodesCoord] = yb;
        NrOfNewNodesCoord ++;
    }
//     printf("NrOfNewNodesCoord = %ld\n",NrOfNewNodesCoord); 

    *NrOfNewNodesCoord_ptr = NrOfNewNodesCoord;
    *NewNodesCoordX_ptr = NewNodesCoordX;
    *NewNodesCoordY_ptr = NewNodesCoordY;
    
    return SUCCESS;
}
