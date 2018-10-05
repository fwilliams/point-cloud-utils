/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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

#ifndef __VCG_GL_MESH_ATTRIBUTES_INFO
#define __VCG_GL_MESH_ATTRIBUTES_INFO

#include <vector>
#include <string>
#include <bitset>

namespace vcg
{
	struct GLMeshAttributesInfo
    {
        struct Exception : public std::exception
        {
            Exception(const char* text)
                :std::exception(),_text(text) {}

            ~Exception() throw() {}
            inline const char* what() const throw() {return _text.c_str();}
        private:
            std::string _text;
        };


        /*WARNING!!!! why not a plain and simple enum? because i need to add further values to this enumeration, but the user of the class should not be interested and/or directly manage those other additional values*/
        /*the struct is extended in GLMeshAttributesFeeder class, introducing the ATT_VERTINDICES and ATT_EDGEINDICES values for the class internal use*/ 
        struct ATT_NAMES
        {
            static const unsigned int ATT_VERTPOSITION = 0;
            static const unsigned int ATT_VERTNORMAL = 1;
            static const unsigned int ATT_FACENORMAL = 2;
            static const unsigned int ATT_VERTCOLOR = 3;
            static const unsigned int ATT_FACECOLOR = 4;
            static const unsigned int ATT_VERTTEXTURE = 5;
            static const unsigned int ATT_WEDGETEXTURE = 6;
            enum {ATT_ARITY = 7};


            ATT_NAMES()
                :_val(ATT_VERTPOSITION)
            {
            }

            ATT_NAMES(unsigned int att)
            {
                if (att >= ATT_NAMES::enumArity())
                    throw Exception("Out of range value\n");
                else
                    _val = att;
            }

            static unsigned int enumArity()
            {
                return ATT_NAMES::ATT_ARITY;
            }

            operator unsigned int() const
            {
                return _val;
            }

            bool operator==(unsigned int r) const
            {
                return (_val == r);
            }

            bool operator!=(unsigned int r) const
            {
                return (_val != r);
            }

        protected:
            unsigned int _val;
        };

        enum PRIMITIVE_MODALITY
        {            
            PR_POINTS = 0,              
            PR_WIREFRAME_EDGES = 1,      
            PR_WIREFRAME_TRIANGLES = 2, 
            PR_SOLID = 3,               
            PR_ARITY = 4                 
        };

        static PRIMITIVE_MODALITY next(PRIMITIVE_MODALITY pm)
        {
            int tmp = static_cast<int>(pm);
            if (tmp == PR_ARITY)
                throw Exception("PRIMITIVE_MODALITY iterator: PR_ARITY passed as parameter!");
            ++tmp;          
            return static_cast<PRIMITIVE_MODALITY>(tmp);
        }

        typedef std::bitset<PR_ARITY> PRIMITIVE_MODALITY_MASK;

        template<typename ATT_NAMES_DERIVED_CLASS>
        class RenderingAtts
        {
        public:
            RenderingAtts(bool defaultvalue = false)
            {
                reset(defaultvalue);
            }

            RenderingAtts(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& att)
            {
                reset();
                //_atts = new bool[ATT_NAMES_DERIVED_CLASS::enumArity()];
                for(unsigned int ii = 0;ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                    (*this)[ii] = att[ii];
            }

            ~RenderingAtts()
            {
            }

            RenderingAtts<ATT_NAMES_DERIVED_CLASS>& operator=(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& att)
            {
                reset();
                for(unsigned int ii = 0;ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                    (*this)[ii] = att[ii];
                return (*this);
            }

            /*bool operator[](ATT_NAMES_DERIVED_CLASS att) const
            {
                unsigned int ii = att;
                if (ii >= ATT_NAMES_DERIVED_CLASS::enumArity())
                    throw GLFeederException("Out of range value\n");
                return _atts[ii];
            }

            bool& operator[](ATT_NAMES_DERIVED_CLASS att)
            {
                unsigned int ii = att;
                if (ii >= ATT_NAMES_DERIVED_CLASS::enumArity())
                    throw GLFeederException("Out of range value\n");
                return _atts[ii];
            }*/

            bool operator[](unsigned int ind) const
            {
                if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
                    throw Exception("Out of range value\n");
                return _atts[ind];
            }

            bool& operator[](unsigned int ind)
            {
                if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
                    throw Exception("Out of range value\n");
                return _atts[ind];
            }

            void reset(bool defaultvalue = false)
            {
                //delete[] _atts;
                //_atts = new bool[ATT_NAMES_DERIVED_CLASS::enumArity()];
                for(unsigned int ii = 0;ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                    _atts[ii] = defaultvalue;
            }

            static RenderingAtts<ATT_NAMES_DERIVED_CLASS> unionSet(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
            {
                RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
                for(unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                    res[ii] = a[ii] || b[ii];
                return res;
            }

            static RenderingAtts<ATT_NAMES_DERIVED_CLASS> complementSet(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
            {
                /*TRUTH TABLE*/
                //this[ATT_NAMES] | rq[ATT_NAMES] | res
                //    true        |     true      | false
                //    true        |     false     | true
                //    false       |     true      | false
                //    false       |     false     | false

                RenderingAtts<ATT_NAMES_DERIVED_CLASS> res = a;
                for(unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                {
                    if (res[ii])
                        res[ii] = !(b[ii]);
                }

                return res;
            }

            static RenderingAtts<ATT_NAMES_DERIVED_CLASS> intersectionSet(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
            {
                RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
                for(unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity();++ii)
                    res[ii] = a[ii] && b[ii];
                return res;
            }

            //template<typename MESHTYPE>
            //static void computeARequestedAttributesSetCompatibleWithMesh(const MESHTYPE& mesh,const PRIMITIVE_MODALITY_MASK,RenderingAtts<ATT_NAMES_DERIVED_CLASS>& rqatt)
            //{
            //    if (mesh.VN() == 0)
            //    {
            //        rqatt.reset();
            //        return;
            //    }

            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTPOSITION] = true;
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTNORMAL] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(mesh);
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FACENORMAL] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FACENORMAL] && vcg::tri::HasPerFaceNormal(mesh);
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTCOLOR] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(mesh);
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FACECOLOR] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FACECOLOR] && vcg::tri::HasPerFaceColor(mesh);
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FIXEDCOLOR] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_FIXEDCOLOR];
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTTEXTURE] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_VERTTEXTURE] && vcg::tri::HasPerVertexTexCoord(mesh);
            //    rqatt[ATT_NAMES_DERIVED_CLASS::ATT_WEDGETEXTURE] = rqatt[ATT_NAMES_DERIVED_CLASS::ATT_WEDGETEXTURE] && vcg::tri::HasPerWedgeTexCoord(mesh);
            //}
        protected:
            /*an array of enumArity() bool values*/
            bool _atts[ATT_NAMES_DERIVED_CLASS::ATT_ARITY];
            //std::bitset<ATT_NAMES_DERIVED_CLASS::ATT_ARITY> _atts;
        };

        typedef RenderingAtts<ATT_NAMES> RendAtts;

        struct DebugInfo
        {
            std::string _tobeallocated;
            std::string _tobedeallocated;
            std::string _tobeupdated;

            std::string _currentlyallocated;

            std::vector<std::string> _perviewdata;

            DebugInfo()
                :_tobeallocated(),_tobedeallocated(),_tobeupdated(),_currentlyallocated(),_perviewdata()
            {

            }

            void reset()
            {
                _tobeallocated.clear();
                _tobedeallocated.clear();
                _tobeupdated.clear();
                _currentlyallocated.clear();
                _perviewdata.clear();
            }

            static const char* primitiveName(size_t ind)
            {
                static std::string res;

                if (ind == size_t(PR_POINTS))
                    res = std::string("PR_POINTS");

                if (ind == size_t(PR_WIREFRAME_EDGES))
                    res = std::string("PR_WIREFRAME_EDGES");

                if (ind == size_t(PR_WIREFRAME_TRIANGLES))
                    res = std::string("PR_WIREFRAME_TRIANGLES");

                if (ind == size_t(PR_SOLID))
                    res = std::string("PR_SOLID");

                return res.c_str();
            }
        };


    protected:
        struct INT_ATT_NAMES : public ATT_NAMES
        {
            /*WARNING!!!!!! the edges index bo it's just used only by the edges and quads meshes, NOT by the triangle meshes. Triangles meshes use just the vertex index array*/
            /*WHY? cause quads meshes need both index arrays. One to render the "usual" mesh, the other one to render the wireframe quadrangulation on top of it*/
            /*A triangles meshes rendered in wireframe or in solid wireframe, in order to save GPU memory, use the glPolygonMode approach*/
            /*Edges meshes uses the edges index array just for a matter of coherence*/ 

            /*WARNING!!!! to be changed whit ATT_NAMES::enumArity() (and so on...) as soon as constexpr will be supported by most of the old c++ compilers*/
            static const unsigned int ATT_VERTINDICES = 7;
            static const unsigned int ATT_EDGEINDICES = 8;
            enum {ATT_ARITY = 9};

            INT_ATT_NAMES()
                :ATT_NAMES()
            {
            }

            INT_ATT_NAMES(unsigned int att)
                :ATT_NAMES()
            {
                if (att >= INT_ATT_NAMES::enumArity())
                    throw Exception("Out of range value\n");
                else
                    _val = att;
            }

            static unsigned int enumArity() 
            {
                return INT_ATT_NAMES::ATT_ARITY;
            }

            operator unsigned int() const
            {
                return _val;
            }
        };

        class InternalRendAtts : public RenderingAtts<INT_ATT_NAMES>
        {
        public:
            InternalRendAtts()
                :RenderingAtts<INT_ATT_NAMES>()
            {
            }

            InternalRendAtts(const RendAtts& reqatt)
                :RenderingAtts<INT_ATT_NAMES>()
            {
                
                for(unsigned int ii = 0;ii < ATT_NAMES::enumArity();++ii)
                {
                    (*this)[ii] = reqatt[ii];
                }

                (*this)[INT_ATT_NAMES::ATT_VERTINDICES] = false;
                (*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = false;
            }

            InternalRendAtts(const RendAtts& reqatt,PRIMITIVE_MODALITY pm)
                :RenderingAtts<INT_ATT_NAMES>()
            {
                for(unsigned int ii = 0;ii < ATT_NAMES::enumArity();++ii)
                    (*this)[ii] = reqatt[ii];

                (*this)[INT_ATT_NAMES::ATT_VERTINDICES] = isVertexIndexingRequired(reqatt,pm);
                (*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = isEdgeIndexingRequired(pm);
            }


            InternalRendAtts(const RenderingAtts<INT_ATT_NAMES>& r)
            :RenderingAtts<INT_ATT_NAMES>(r)
            {
            }

            //upcast from InternalRendAtts to RendAtts
            operator RendAtts() const
            {
                RendAtts rendatt;
                for(unsigned int ii = 0;ii < ATT_NAMES::enumArity();++ii)
                    rendatt[ii] = _atts[ii];
                return rendatt;
            }

            InternalRendAtts& setIndexingIfNeeded(PRIMITIVE_MODALITY pm)
            {
                (*this)[INT_ATT_NAMES::ATT_VERTINDICES] = isVertexIndexingRequired((*this),pm);
                (*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = isEdgeIndexingRequired(pm);
                return (*this);
            }

            static bool isPerVertexAttribute(INT_ATT_NAMES name)
            {
                return ((name == INT_ATT_NAMES::ATT_VERTPOSITION) || (name == INT_ATT_NAMES::ATT_VERTNORMAL) || (name == INT_ATT_NAMES::ATT_VERTCOLOR) || (name == INT_ATT_NAMES::ATT_VERTTEXTURE));
            }

            static bool replicatedPipelineNeeded(const RendAtts& rqatt)
            {
                return (rqatt[INT_ATT_NAMES::ATT_FACENORMAL] || rqatt[INT_ATT_NAMES::ATT_FACECOLOR] || rqatt[INT_ATT_NAMES::ATT_WEDGETEXTURE]);
            }

            static bool isVertexIndexingRequired(const RendAtts& rqatt,PRIMITIVE_MODALITY pm)
            {
                return (!replicatedPipelineNeeded(rqatt) && ((pm == PR_SOLID) || (pm == PR_WIREFRAME_TRIANGLES)));
            }

            static bool isEdgeIndexingRequired(PRIMITIVE_MODALITY pm)
            {
                return (pm == PR_WIREFRAME_EDGES);
            }

            /*static void suggestedMinimalAttributeSetForPrimitiveModalityMask(PRIMITIVE_MODALITY_MASK pm,RenderingAtts<INT_ATT_NAMES>& atts)
            {
                if ((pm == (unsigned int)(PR_NONE)) || (pm == (unsigned int)(PR_BBOX)))
                {
                    atts.reset();
                    return;
                }

                if (pm & PR_POINTS)
                {
                    atts[INT_ATT_NAMES::ATT_VERTPOSITION] = true;
                }

                if (pm & PR_WIREFRAME_EDGES)
                {
                    atts[INT_ATT_NAMES::ATT_VERTPOSITION] = true;
                    atts[INT_ATT_NAMES::ATT_EDGEINDICES] = true;
                }

                if ((pm & PR_WIREFRAME_TRIANGLES) || (pm & PR_SOLID))
                {
                    atts[INT_ATT_NAMES::ATT_VERTPOSITION] = true;
                    atts[INT_ATT_NAMES::ATT_VERTINDICES] = true;
                }

            }*/
        };
    };
}

#endif
