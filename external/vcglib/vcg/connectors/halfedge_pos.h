/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef VCG_HEDGE_POS
#define VCG_HEDGE_POS

namespace vcg
{
    namespace hedge
    {
        /*!
          * \brief Class implementing a Pos using halfedges
          *
          */
        template <class MeshType> class Pos
        {

        public:

            typedef typename MeshType::VertexPointer VertexPointer;
            typedef typename MeshType::EdgePointer EdgePointer;
            typedef typename MeshType::HEdgePointer HEdgePointer;
            typedef typename MeshType::FacePointer FacePointer;

            /*!
             * Halfedge used to move in the mesh
             */
            HEdgePointer he;

            /*!
             * Direction of movement:
             *
             * 0 = clockwise
             *
             * 1 = counter-clockwise
             */
            bool direction;

            /*!
             *
             */
            Pos(HEdgePointer hep, bool dir)
            {
                he = hep;
                direction = dir;
            }

            /*!
             *
             */
            Pos(HEdgePointer hep)
            {
                he = hep;
                direction = 1;
            }

            /*!
             * Changes vertex mantaininge the same edge and the same face
             */
            void FlipV()
            {
                direction = !direction;
            }

            /*!
             * Changes edge and hedge mantaining the same vertex and the same face
             */
            void FlipE()
            {
                if(!direction)
                    he = he->HNp();
                else
                    if(he->HasHPrevAdjacency())
                        he = he->HPp();
                    else
                    {
                        HEdgePointer aux = he;
                        while(aux->HNp() != he)
                            aux = aux->HNp();
                        he = aux;
                    }

                direction = !direction;
            }

            /*!
             * Changes face mantaininge the same vertex and the same edge
             */
            void FlipF()
            {
                direction = !direction;
                he = he->HOp();
            }

            /*!
             * Gets pointed vertex
             */
            VertexPointer V()
            {
                if(direction)
                    return he->HVp();
                else
                    return he->HOp()->HVp();
            }

			/*!
             * Gets opposite vertex
             */
            VertexPointer Vo()
            {
                if(!direction)
                    return he->HVp();
                else
                    return he->HOp()->HVp();
            }
			
            /*!
             * Gets pointed hedge
             */
            HEdgePointer HE()
            {
                return he;
            }

            /*!
             * Gets pointed edge
             */
            EdgePointer E()
            {
                return he->HEp();
            }

            /*!
             * Gets pointed face
             */
            FacePointer F()
            {
                return he->HFp();
            }

            /*!
             * Operator to check if two Pos are equal
             */
            inline bool operator == ( Pos const & p ) const
            {
                return (he == p.he && direction == p.direction);
            }

            /*!
             * Operator to check if two Pos are different
             */
            inline bool operator != ( Pos const & p ) const
            {
                return (he != p.he || direction != p.direction);
            }

        };

    }
}

#endif // VCG_HEDGE_POS

