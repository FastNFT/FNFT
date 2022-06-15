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

/**
 * @file fnft__delaunay_triangulation.h
 * @brief Delaunay triangulation.
 * @ingroup delaunay_triangulation
 */

#ifndef FNFT__DELAUNAY_TRIANGULATION_H
#define FNFT__DELAUNAY_TRIANGULATION_H

#include <stdio.h>
#include "fnft__misc.h"
#include <string.h> // for memcpy
#include "fnft__errwarn.h"

/**
 * @struct fnft__delaunay_triangulation_data_t;
 * @brief Stores the data from Delaunay triangulation. The structures are 
 * passed between the routines \link fnft_nsev \endlink and \link fnft__delaunay_triangulation \endlink.
 * @ingroup fnft
 * @ingroup data_types
 *
 * @var fnft__delaunay_triangulation_data_t::NodesMax
 *  This is the maximum number of nodes that will be triangulated. The
 *  memory is allocated based on this value.
 * 
 * @var fnft__delaunay_triangulation_data_t::NodeOfTriangles1
 *  This is the pointer to the first node/vertex of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::NodeOfTriangles2
 *  This is the pointer to the second node/vertex of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::NodeOfTriangles3
 *  This is the pointer to the third node/vertex of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::EdgesOfTriangles1
 *  This is the pointer to the first edge of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::EdgesOfTriangles2
 *  This is the pointer to the second edge of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::EdgesOfTriangles3
 *  This is the pointer to the third edge of the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::Edges1
 *  This is the pointer to the first node of the edges.
 *
 * @var fnft__delaunay_triangulation_data_t::Edges2
 *  This is the pointer to the second node of the edges.
 *  NOTE: Edges1[i]<Edges2[i].
 *
 * @var fnft__delaunay_triangulation_data_t::NodesCoordX
 * This is the pointer to the x-coordinates of all the triangulated nodes.
 *
 * @var fnft__delaunay_triangulation_data_t::NodesCoordY
 * This is the pointer to the y-coordinates of all the triangulated nodes.
 *
 * @var fnft__delaunay_triangulation_data_t::NrOfNodes
 * This is the number of currently triangulated nodes.
 *
 * @var fnft__delaunay_triangulation_data_t::NrOfTriangles
 * This is the number of currently built triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::NrOfNodesToRemove
 * This is the number of extra nodes used for initialization of the 
 * triangluation.
 *
 * @var fnft__delaunay_triangulation_data_t::NrOfEdges
 * This is the number of currently built edges.
 *
 * @var fnft__delaunay_triangulation_data_t::XCoordOfCCOfTriangles
 * This is the pointer to the x-coordinates of the circumcenters of 
 * all the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::YCoordOfCCOfTriangles
 * This is the pointer to the y-coordinates of the circumcenters of 
 * all the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::RadiusOfCCOfTriangles
 * This is the pointer to the radii of the circumcircles of 
 * all the triangles.
 *
 * @var fnft__delaunay_triangulation_data_t::StatusOfTriangles
 * This is flag array which indicates whether a triangle i is valid if 
 * StatusOfTriangles[i] == 1.
 *
 * @var fnft__delaunay_triangulation_data_t::CheckEdge
 * This is flag array which indicates whether a for a triangle its edges
 * need to be checked against list of edges and added if not in there already.
 *
 * @var fnft__delaunay_triangulation_data_t::RemoveSuperRectangle
 * This is a flag used during triangluation.
 * If RemoveSuperRectangle=1 then the triangles which have the initialization
 * nodes are set as invalid.
 *
 */
typedef struct {
    FNFT_UINT NodesMax;
    FNFT_UINT * NodeOfTriangles1;
    FNFT_UINT * NodeOfTriangles2;
    FNFT_UINT * NodeOfTriangles3;
    FNFT_UINT * EdgesOfTriangles1;
    FNFT_UINT * EdgesOfTriangles2;
    FNFT_UINT * EdgesOfTriangles3;
    FNFT_UINT * Edges1;
    FNFT_UINT * Edges2;
    FNFT_REAL * NodesCoordX;
    FNFT_REAL * NodesCoordY;
    FNFT_UINT NrOfNodes;
    FNFT_UINT NrOfTriangles;
    FNFT_UINT NrOfNodesToRemove;
    FNFT_UINT NrOfEdges;
    FNFT_REAL * XCoordOfCCOfTriangles;
    FNFT_REAL * YCoordOfCCOfTriangles;
    FNFT_REAL * RadiusOfCCOfTriangles;
    FNFT_UINT * StatusOfTriangles;
    FNFT_UINT * CheckEdge;
    FNFT_UINT RemoveSuperRectangle;
}fnft__delaunay_triangulation_data_t;


/**
 * @brief Peforms Delaunay triangulation.
 * 
 * The function performs Delaunay triangulation and is based on the iterative
 * Bowyer–Watson algorithm. The iterative nature allows for addition of new points
 * even in separate calls to the routine. The actual implementation is based on
 * the reference
 *      - Sloan and Houlsby, <a href="https://doi.org/10.1016/0141-1195(84)90003-2">&quot; An implementation of Watson's algorithm for computing 2-dimensional delaunay triangulations,&quot;</a> Advances in Engineering Software, 6(4), pp. 192--197, 1978.
 *
 * @param[in] DT_ptr It the pointer to a structure of type  \link fnft__delaunay_triangulation_data_t \endlink.
 * @param[in] NodesMax The maximum number of nodes that will be triangulated.
 * @param[in] NrOfNewNodes The new number of nodes that are to be added to the triangulation.
 * @param[in] NewNodesCoordX Array of length NrOfNewNodes with the x-coordinates of the new
 * points that are to be added to the triangulation.
 * @param[in] NewNodesCoordY Array of length NrOfNewNodes with the y-coordinates of the new
 * points that are to be added to the triangulation.
 * @param[in, out] status_flag_ptr Pointer to the status_flag. Should be 0 on entry.
 * If a duplicate point is being added status_flag is set to 1 and the function 
 * stops triangulation and returns \link FNFT_SUCCESS \endlink. Similarly
 * status_flag is set to 2 if the routine runs out of memory for triangles
 * and set to 3 if the routine runs out of memory for edges.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup delaunay_triangulation
 */
FNFT_INT fnft__delaunay_triangulation(fnft__delaunay_triangulation_data_t * DT_ptr,
        FNFT_UINT NodesMax, FNFT_UINT NrOfNewNodes, FNFT_REAL * NewNodesCoordX, 
        FNFT_REAL * NewNodesCoordY, FNFT_UINT * status_flag_ptr);

/**
 * @brief Finds all the triangles which share specified edges.
 * 
 * The function finds all the triangles which have the listed edges.
 * The triangles are listed in the form of their three vertices. The edges
 * are listed as their two end points.
 *
 * @param[in] NrOfElements Number of triangle elements.
 * @param[in] Elements1 Array of length NrOfElements which contain the first
 * vertex of the triangles.
 * @param[in] Elements2 Array of length NrOfElements which contain the second
 * vertex of the triangles.
 * @param[in] Elements3 Array of length NrOfElements which contain the third
 * vertex of the triangles.
 * @param[in] NrOfNodes Number of edges.
 * @param[in] Nodes1 Array of length NrOfNodes which contain the first
 * end point of the edges.
 * @param[in] Nodes2 Array of length NrOfNodes which contain the second
 * end point of the edges.
 * @param[out] NrOfCandidateElements_ptr This should be the pointer to  NrOfCandidateElements
 * which on return will contain the length of ArrayOfCandidateElements.
 * @param[out] ArrayOfCandidateElements_ptr This is the pointer to the ArrayOfCandidateElements array
 * which contains the indices of triangles which were found to have the specified edges. 
 * NOTE: The indices of the triangles may have repitition if the triangle has
 * more than one of the specified edges. Atleast 3*NrOfElements amount of
 * memory should be alloted.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup delaunay_triangulation
 */
FNFT_INT fnft__delaunay_triangulation_edgeAttachments(FNFT_UINT NrOfElements,
        FNFT_UINT *Elements1, FNFT_UINT *Elements2, FNFT_UINT *Elements3,
        FNFT_UINT NrOfNodes, FNFT_UINT *Nodes1, FNFT_UINT *Nodes2,
        FNFT_UINT *NrOfCandidateElements_ptr,
        FNFT_UINT *ArrayOfCandidateElements);

/**
 * @brief Finds all the triangles which have specified vertices.
 * 
 * The function finds all the triangles which have the listed vertices.
 * The triangles are listed in the form of their three vertices. 
 *
 * @param[in] NrOfElements Number of triangle elements.
 * @param[in] Elements1 Array of length NrOfElements which contain the first
 * vertex of the triangles.
 * @param[in] Elements2 Array of length NrOfElements which contain the second
 * vertex of the triangles.
 * @param[in] Elements3 Array of length NrOfElements which contain the third
 * vertex of the triangles.
 * @param[in] NrOfNodes Number of nodes/vertices.
 * @param[in] Nodes Array of length NrOfNodes which contains the nodes.
 * @param[out] NrOfCandidateElements_ptr This should be the pointer to  NrOfCandidateElements
 * which on return will contain the length of ArrayOfCandidateElements.
 * @param[out] ArrayOfCandidateElements_ptr This is the pointer to the ArrayOfCandidateElements array
 * which contains the indices of triangles which were found to have the specified nodes.
 * NOTE: The indices of the triangles may have repitition if the triangle has
 * more than one of the specified nodes. Atleast 3*NrOfElements amount of
 * memory should be alloted.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup delaunay_triangulation
 */
FNFT_INT fnft__delaunay_triangulation_vertexAttachments(FNFT_UINT NrOfElements,
        FNFT_UINT *Elements1, FNFT_UINT *Elements2, FNFT_UINT *Elements3,
        FNFT_UINT NrOfNodes, FNFT_UINT *Nodes,
        FNFT_UINT *NrOfCandidateElements_ptr,
        FNFT_UINT *ArrayOfCandidateElements);

/**
 * @brief Generates rectangular mesh.
 * 
 * The function generates a mesh for a rectangular domain based on limits of 
 * the domain and chosen step-size.
 * @param[in] xb Real scalar value which is the smallest value of the x-coordinate.
 * @param[in] xe Real scalar value which is the largest value of the x-coordinate.
 * @param[in] yb Real scalar value which is the smallest value of the y-coordinate.
 * @param[in] ye Real scalar value which is the largest value of the y-coordinate.
 * @param[in] r Real scalar step size of the mesh to be generated. 
 * @param[out] NrOfNewNodesCoord_ptr This should be the pointer to NrOfNewNodesCoord
 * which on return will contain the number of points in the generated mesh.
 * @param[out] NewNodesCoordX_ptr This is the pointer to the NewNodesCoordX array
 * which contains the x-coordinates of the generated mesh points. The memory for the
 * array is allocated by the routine and should be freed by the calling function.
 * @param[out] NewNodesCoordY_ptr This is the pointer to the NewNodesCoordY array
 * which contains the y-coordinates of the generated mesh points. The memory for the
 * array is allocated by the routine and should be freed by the calling function.
 * @return \link FNFT_SUCCESS \endlink or one of the FNFT_EC_... error codes
 *  defined in \link fnft_errwarn.h \endlink.
 * @ingroup delaunay_triangulation
 */
FNFT_INT fnft__delaunay_triangulation_rect_dom(FNFT_REAL xb, FNFT_REAL xe, FNFT_REAL yb, FNFT_REAL ye, FNFT_REAL r,
        FNFT_UINT * const NrOfNewNodesCoord_ptr, FNFT_REAL ** NewNodesCoordX_ptr, FNFT_REAL ** NewNodesCoordY_ptr);
        
        
#ifdef FNFT_ENABLE_SHORT_NAMES
#define  delaunay_triangulation(...) fnft__delaunay_triangulation(__VA_ARGS__)
#define  delaunay_triangulation_vertexAttachments(...) fnft__delaunay_triangulation_vertexAttachments(__VA_ARGS__)
#define  delaunay_triangulation_edgeAttachments(...) fnft__delaunay_triangulation_edgeAttachments(__VA_ARGS__)
#define  delaunay_triangulation_rect_dom(...) fnft__delaunay_triangulation_rect_dom(__VA_ARGS__)
#endif

#endif
