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

#ifndef __VCG_GL_MESH_ATTRIBUTES_FEEDER
#define __VCG_GL_MESH_ATTRIBUTES_FEEDER

#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <climits>
#include <string>

#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/space/color4.h>
#include <vcg/math/matrix44.h>
#include<wrap/system/memory_info.h>
#include <wrap/gl/gl_mesh_attributes_info.h>


namespace vcg
{
    struct RenderingModalityGLOptions
    {
        bool _perbbox_enabled;

        bool _perbbox_fixed_color_enabled;
        bool _perpoint_fixed_color_enabled;
        bool _perwire_fixed_color_enabled;
        bool _persolid_fixed_color_enabled;

        Color4b _perbbox_fixed_color;
        Color4b _perpoint_fixed_color;
        Color4b _perwire_fixed_color;
        Color4b _persolid_fixed_color;

        bool _perbbox_mesh_color_enabled;
        bool _perpoint_mesh_color_enabled;
        bool _perwire_mesh_color_enabled;
        bool _persolid_mesh_color_enabled;
        
        bool _perpoint_noshading;
        bool _perwire_noshading;
        bool _persolid_noshading;

        bool _perpoint_dot_enabled;

        float _perpoint_pointsize;
        bool _perpoint_pointsmooth_enabled;
        bool _perpoint_pointattenuation_enabled;

        float _perwire_wirewidth;

        RenderingModalityGLOptions()
        {
            _perbbox_enabled = false;

            _perbbox_fixed_color_enabled = true;
            _perpoint_fixed_color_enabled = false;
            _perwire_fixed_color_enabled = true;
			_persolid_fixed_color_enabled = true;

            _perbbox_fixed_color = vcg::Color4b(Color4b::White);
            _perpoint_fixed_color = vcg::Color4b(Color4b::White);
            _perwire_fixed_color = Color4b(Color4b::DarkGray);
            _persolid_fixed_color = vcg::Color4b(Color4b::White);

            _perbbox_mesh_color_enabled = false;
            _perpoint_mesh_color_enabled = false;
            _perwire_mesh_color_enabled = false;
			_persolid_mesh_color_enabled = false;

            _perpoint_dot_enabled = false;

            _perpoint_noshading = false;
            _perwire_noshading = true;
            _persolid_noshading = false;

            _perpoint_pointsize = 1.0f;
            _perpoint_pointsmooth_enabled = false;
            _perpoint_pointattenuation_enabled = true;
            _perwire_wirewidth = 1.0f;
        }

        RenderingModalityGLOptions(const RenderingModalityGLOptions& opts)
        {
            copyData(opts);
        }

        virtual ~RenderingModalityGLOptions()
        {
        }

        RenderingModalityGLOptions& operator=(const RenderingModalityGLOptions& opts)
        {
            copyData(opts);
            return (*this);
        }


    private:
        void copyData(const RenderingModalityGLOptions& opts)
        {
            _perbbox_enabled = opts._perbbox_enabled;

            _perpoint_dot_enabled = opts._perpoint_dot_enabled; 
            _perpoint_pointsize = opts._perpoint_pointsize;
            _perpoint_pointsmooth_enabled = opts._perpoint_pointsmooth_enabled;
            _perpoint_pointattenuation_enabled = opts._perpoint_pointattenuation_enabled;

            _perbbox_fixed_color_enabled = opts._perbbox_fixed_color_enabled;
            _perpoint_fixed_color_enabled = opts._perpoint_fixed_color_enabled;
            _perwire_fixed_color_enabled = opts._perwire_fixed_color_enabled;
            _persolid_fixed_color_enabled = opts._persolid_fixed_color_enabled;

            _perbbox_mesh_color_enabled = opts._perbbox_mesh_color_enabled;
            _perpoint_mesh_color_enabled = opts._perpoint_mesh_color_enabled;
            _perwire_mesh_color_enabled = opts._perwire_mesh_color_enabled;
            _persolid_mesh_color_enabled = opts._persolid_mesh_color_enabled;
            
            _perbbox_fixed_color = opts._perbbox_fixed_color;
            _perpoint_fixed_color = opts._perpoint_fixed_color;
            _perwire_fixed_color = opts._perwire_fixed_color;
            _persolid_fixed_color = opts._persolid_fixed_color;

            _perpoint_noshading = opts._perpoint_noshading;
            _perwire_noshading = opts._perwire_noshading;
            _persolid_noshading = opts._persolid_noshading;

            _perwire_wirewidth = opts._perwire_wirewidth;
        }
    };

    template<typename GL_OPTIONS_DERIVED_TYPE = RenderingModalityGLOptions>
    class PerViewData : public GLMeshAttributesInfo
    {
    public:
        PerViewData()
            :_pmmask(),_intatts(PR_ARITY),_glopts(NULL)
        {
            reset();
        }


        PerViewData(const PerViewData<GL_OPTIONS_DERIVED_TYPE>& dt)
            :_pmmask(dt._pmmask),_intatts(dt._intatts),_glopts(NULL)
        {
            if (dt._glopts != NULL)
                _glopts = new GL_OPTIONS_DERIVED_TYPE(*(dt._glopts));
        }

        ~PerViewData()
        {
            _intatts.clear();
            delete _glopts;
        }

        PerViewData& operator=(const PerViewData<GL_OPTIONS_DERIVED_TYPE>& dt)
        {
            _pmmask = dt._pmmask;
            _intatts = dt._intatts;
            if (dt._glopts != NULL)
                _glopts = new GL_OPTIONS_DERIVED_TYPE(*(dt._glopts));
            return (*this);
        }

        bool set(PRIMITIVE_MODALITY pm,const RendAtts& atts)
        {
            size_t pmind(pm);
            if (pm >= _intatts.size())
                return false;
            //_pmmask.set(pm);
            _intatts[pmind] = InternalRendAtts(atts,pm);
            _pmmask.set(size_t(pm),_intatts[pmind][INT_ATT_NAMES::ATT_VERTPOSITION]);
            return true;
        }

        bool set(PRIMITIVE_MODALITY pm,ATT_NAMES att,bool onoff)
        {
            size_t pmind(pm);
            if (pm >= _intatts.size())
                return false;
            _intatts[pmind][att] = onoff;
            _pmmask.set(size_t(pm),_intatts[pmind][INT_ATT_NAMES::ATT_VERTPOSITION]);
            if (_pmmask.test(size_t(pm)))
                _intatts[pmind].setIndexingIfNeeded(pm);
            return true;
        }

        bool set(PRIMITIVE_MODALITY pm,bool onoff)
        {
            return set(pm,INT_ATT_NAMES::ATT_VERTPOSITION,onoff);
        }

        void set(const GL_OPTIONS_DERIVED_TYPE& opts)
        {
            delete _glopts;
            _glopts = new GL_OPTIONS_DERIVED_TYPE(opts); 
        }

        bool isPrimitiveActive(PRIMITIVE_MODALITY pm) const
        {
            if (pm == PR_ARITY)
                return false;
            return (_pmmask.test(pm) && _intatts[size_t(pm)][INT_ATT_NAMES::ATT_VERTPOSITION]);
        }

        PRIMITIVE_MODALITY_MASK getPrimitiveModalityMask() const
        {
            return _pmmask;
        }

        bool get(PRIMITIVE_MODALITY pm,RendAtts& atts) const
        {
            size_t pmind(pm);
            if (pm >= _intatts.size())
                return false;
            atts = _intatts[pmind];
            return true;
        } 

        bool get(GL_OPTIONS_DERIVED_TYPE& opts) const
        {
            if (_glopts == NULL)
                return false;
            opts = (*_glopts);
            return true;
        }

        void reset(bool deleteglopts = true)
        {
            _pmmask.reset();
            for(typename PerRendModData::iterator it = _intatts.begin();it != _intatts.end();++it)
                it->reset();
            if (deleteglopts)
            {
                delete _glopts;
                _glopts = 0;
            }
        }

    protected:
        template<typename MESH_TYPE,typename UNIQUE_VIEW_ID_TYPE, typename XX_GL_OPTIONS_DERIVED_TYPE> friend class NotThreadSafeGLMeshAttributesMultiViewerBOManager;

        typedef std::vector<InternalRendAtts> PerRendModData;

        PRIMITIVE_MODALITY_MASK _pmmask;
        PerRendModData _intatts;

        GL_OPTIONS_DERIVED_TYPE* _glopts;
    };

    /****************************************************WARNING!!!!!!!!!!!!!!!!!*********************************************************************************************/ 
    //You must inherit from NotThreadSafeGLMeshAttributesMultiViewerBOManager, providing thread safe mechanisms, in order to use the bo facilities exposed by the class.
    //In wrap/qt/qt_thread_safe_memory_rendering.h you will find a ready to use class based on QtConcurrency module.
    /*************************************************************************************************************************************************************************/
    template<typename MESH_TYPE,typename UNIQUE_VIEW_ID_TYPE = unsigned int,typename GL_OPTIONS_DERIVED_TYPE = RenderingModalityGLOptions>
    class NotThreadSafeGLMeshAttributesMultiViewerBOManager : public GLMeshAttributesInfo
    {
    public:
        typedef PerViewData<GL_OPTIONS_DERIVED_TYPE> PVData;

    protected:
        /****************************************************WARNING!!!!!!!!!!!!!!!!!*********************************************************************************************/ 
        //You must inherit from NotThreadSafeGLMeshAttributesMultiViewerBOManager, providing thread safe mechanisms, in order to use the bo facilities exposed by the class.
        //In wrap/qt/qt_thread_safe_memory_rendering.h you will find a ready to use class based on QtConcurrency module.
        /*************************************************************************************************************************************************************************/

        NotThreadSafeGLMeshAttributesMultiViewerBOManager(/*const*/ MESH_TYPE& mesh,MemoryInfo& meminfo, size_t perbatchprimitives)
            :_mesh(mesh),_gpumeminfo(meminfo),_bo(INT_ATT_NAMES::enumArity(),NULL),_currallocatedboatt(),_perbatchprim(perbatchprimitives),_chunkmap(),_borendering(false),_edge(),_meshverticeswhenedgeindiceswerecomputed(0),_meshtriangleswhenedgeindiceswerecomputed(0),_tr(),_debugmode(false),_loginfo(),_meaningfulattsperprimitive(PR_ARITY,InternalRendAtts())
        {
            _tr.SetIdentity();
            _bo[INT_ATT_NAMES::ATT_VERTPOSITION] = new GLBufferObject(3,GL_FLOAT,GL_VERTEX_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_VERTNORMAL] = new GLBufferObject(3,GL_FLOAT,GL_NORMAL_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_FACENORMAL] = new GLBufferObject(3,GL_FLOAT,GL_NORMAL_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_VERTCOLOR] = new GLBufferObject(4,GL_UNSIGNED_BYTE,GL_COLOR_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_FACECOLOR] = new GLBufferObject(4,GL_UNSIGNED_BYTE,GL_COLOR_ARRAY,GL_ARRAY_BUFFER);
            /*MESHCOLOR has not a buffer object associated with it. It's just a call to glColor3f. it's anyway added to the _bo arrays for sake of coherence*/
            //_bo[INT_ATT_NAMES::ATT_FIXEDCOLOR] = NULL;
            _bo[INT_ATT_NAMES::ATT_VERTTEXTURE] = new GLBufferObject(2,GL_FLOAT,GL_TEXTURE_COORD_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_WEDGETEXTURE] = new GLBufferObject(2,GL_FLOAT,GL_TEXTURE_COORD_ARRAY,GL_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_VERTINDICES] = new GLBufferObject(3,GL_UNSIGNED_INT,GL_ELEMENT_ARRAY_BUFFER);
            _bo[INT_ATT_NAMES::ATT_EDGEINDICES] = new GLBufferObject(2,GL_UNSIGNED_INT,GL_ELEMENT_ARRAY_BUFFER);

            initMeaningfulAttsMask();
        }


        ~NotThreadSafeGLMeshAttributesMultiViewerBOManager()
        {
            _edge.clear();
            for(size_t ii = 0;ii < _bo.size();++ii)
                delete _bo[ii];
            _bo.clear();
        }

        /*MeshAttributesUpdate will force the buffer allocation only of the bo rendered at least by one viewer. */
        /*If a filter add to a mesh, for instance, a per vertex color attribute that was not previously rendered, the meshAttributesUpdate() will ignore the attribute until a viewer require explicitly to render the per-vertex-color, too*/
        /*In order to do it, please, call the setPerViewRendAtts() setting up for at least one the existing viewer (or adding a new viewer) in the RendAtts reqatts parameter the reqatts[ATT_VERTCOLOR] to true value. */  
        void meshAttributesUpdated(bool hasmeshconnectivitychanged,const RendAtts& changedrendatts)
        {
            InternalRendAtts tobeupdated(changedrendatts);
            tobeupdated[INT_ATT_NAMES::ATT_VERTINDICES] = hasmeshconnectivitychanged; 
            tobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES] = hasmeshconnectivitychanged;
            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                INT_ATT_NAMES boname(ii);
                if (_bo[boname] != NULL)
                    _bo[boname]->_isvalid = (_bo[boname]->_isvalid) && !(tobeupdated[boname]);    
            }
        }

        bool getPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid,PVData& data) const
        {
            typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
            if (it == _perviewreqatts.end())
                return false;
            data = it->second;
            return true;
        }

        void setPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid,const PVData& data)
        {
            ///cleanup stage...if an attribute impossible for a primitive modality is still here (it should not be...) we change the required atts into the view
            PVData copydt(data);
            for(PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY;pm = next(pm))
                copydt._intatts[pm] = InternalRendAtts::intersectionSet(copydt._intatts[size_t(pm)],_meaningfulattsperprimitive[size_t(pm)]);
            _perviewreqatts[viewid] = copydt;
        }

		void setPerAllViewsInfo(const PVData& data)
		{
			///cleanup stage...if an attribute impossible for a primitive modality is still here (it should not be...) we change the required atts into the view
			PVData copydt(data);
			for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm))
				copydt._intatts[pm] = InternalRendAtts::intersectionSet(copydt._intatts[size_t(pm)], _meaningfulattsperprimitive[size_t(pm)]);
			for (typename ViewsMap::iterator it = _perviewreqatts.begin(); it != _perviewreqatts.end(); ++it)
				it->second = copydt;
		}

        bool removeView(UNIQUE_VIEW_ID_TYPE viewid)
        {
            typename ViewsMap::iterator it = _perviewreqatts.find(viewid);
            if (it == _perviewreqatts.end())
                return false;
            _perviewreqatts.erase(viewid);
            return true;
        }

        void removeAllViews()
        {
            _perviewreqatts.clear();
        }

        void draw(UNIQUE_VIEW_ID_TYPE viewid,const std::vector<GLuint>& textid = std::vector<GLuint>()) const
        {       
            typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
            if (it == _perviewreqatts.end())
                return;

            const PVData& dt = it->second;
            //const InternalRendAtts& atts = it->second._intatts;
			drawFun(dt, textid);
        }

		
		void drawAllocatedAttributesSubset(UNIQUE_VIEW_ID_TYPE viewid,const PVData& dt, const std::vector<GLuint>& textid = std::vector<GLuint>()) const
		{
			typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
			if (it == _perviewreqatts.end())
				return;

			PVData tmp = dt;
			
			if (!(_currallocatedboatt[INT_ATT_NAMES::ATT_VERTPOSITION]))
			{
				for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm))
				{		
					tmp._pmmask[size_t(pm)] = 0;
					tmp._intatts[size_t(pm)] = InternalRendAtts();
				}
			}
			else
			{
				for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm))
				{
					tmp._intatts[size_t(pm)] = InternalRendAtts::intersectionSet(tmp._intatts[size_t(pm)],_meaningfulattsperprimitive[size_t(pm)]);
					tmp._intatts[size_t(pm)] = InternalRendAtts::intersectionSet(tmp._intatts[size_t(pm)],_currallocatedboatt);
				}
			}
			drawFun(dt, textid);
		}

        bool isBORenderingAvailable() const
        {
            return _borendering;
        }

        bool manageBuffers()
        {
            InternalRendAtts tobeallocated;
            InternalRendAtts tobedeallocated;
            InternalRendAtts tobeupdated;
            bool correctlyallocated = false;
            bool arebuffersok = checkBuffersAllocationStatus(tobeallocated,tobedeallocated,tobeupdated);
            if (!arebuffersok)
                correctlyallocated = manageAndFeedBuffersIfNeeded(tobeallocated,tobedeallocated,tobeupdated);
            if (_debugmode)
                debug(tobeallocated,tobedeallocated,tobeupdated);
            return (arebuffersok || correctlyallocated);
        }
 

        void setGLOptions(UNIQUE_VIEW_ID_TYPE viewid,const GL_OPTIONS_DERIVED_TYPE& opts)
        {
            typename ViewsMap::iterator it = _perviewreqatts.find(viewid);
            if (it == _perviewreqatts.end())
                return;
            it->second.set(opts);
        }

        void setTrMatrix(const vcg::Matrix44<typename MESH_TYPE::ScalarType>& tr)
        {
            _tr = tr;
        }

        void setDebugMode(bool isdebug)
        {
            _debugmode = isdebug;
        }

        void getLog(DebugInfo& info)
        {
            info.reset();
            info._tobedeallocated = _loginfo._tobedeallocated;
            info._tobeallocated = _loginfo._tobeallocated;
            info._tobeupdated = _loginfo._tobeupdated;

            info._currentlyallocated = _loginfo._currentlyallocated;
            info._perviewdata = _loginfo._perviewdata;
            _loginfo.reset();
        }

    private:
        void initMeaningfulAttsMask()
        {
            _meaningfulattsperprimitive[PR_POINTS][INT_ATT_NAMES::ATT_VERTPOSITION] = true; 
            _meaningfulattsperprimitive[PR_POINTS][INT_ATT_NAMES::ATT_VERTNORMAL] = true; 
            _meaningfulattsperprimitive[PR_POINTS][INT_ATT_NAMES::ATT_VERTCOLOR] = true; 
            //_meaningfulattsperprimitive[PR_POINTS][INT_ATT_NAMES::ATT_FIXEDCOLOR] = true; 
            _meaningfulattsperprimitive[PR_POINTS][INT_ATT_NAMES::ATT_VERTTEXTURE] = true;

            _meaningfulattsperprimitive[PR_WIREFRAME_EDGES][INT_ATT_NAMES::ATT_VERTPOSITION] = true; 
            _meaningfulattsperprimitive[PR_WIREFRAME_EDGES][INT_ATT_NAMES::ATT_VERTNORMAL] = true; 
            _meaningfulattsperprimitive[PR_WIREFRAME_EDGES][INT_ATT_NAMES::ATT_VERTCOLOR] = true; 
            //_meaningfulattsperprimitive[PR_WIREFRAME_EDGES][INT_ATT_NAMES::ATT_FIXEDCOLOR] = true; 
            _meaningfulattsperprimitive[PR_WIREFRAME_EDGES][INT_ATT_NAMES::ATT_EDGEINDICES] = true;

            _meaningfulattsperprimitive[PR_WIREFRAME_TRIANGLES][INT_ATT_NAMES::ATT_VERTPOSITION] = true; 
            _meaningfulattsperprimitive[PR_WIREFRAME_TRIANGLES][INT_ATT_NAMES::ATT_VERTNORMAL] = true; 
            _meaningfulattsperprimitive[PR_WIREFRAME_TRIANGLES][INT_ATT_NAMES::ATT_VERTCOLOR] = true; 
            //_meaningfulattsperprimitive[PR_WIREFRAME_TRIANGLES][INT_ATT_NAMES::ATT_FIXEDCOLOR] = true;
            _meaningfulattsperprimitive[PR_WIREFRAME_TRIANGLES][INT_ATT_NAMES::ATT_VERTINDICES] = true;

            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_VERTPOSITION] = true; 
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_VERTNORMAL] = true; 
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_FACENORMAL] = true; 
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_VERTCOLOR] = true; 
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_FACECOLOR] = true; 
            //_meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_FIXEDCOLOR] = true; 
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_VERTTEXTURE] = true;
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_WEDGETEXTURE] = true;
            _meaningfulattsperprimitive[PR_SOLID][INT_ATT_NAMES::ATT_VERTINDICES] = true;
        }

        bool hasMeshAttribute(INT_ATT_NAMES attname) const
        {
            switch(attname)
            {
            case(INT_ATT_NAMES::ATT_VERTPOSITION):
                    return true;
            case(INT_ATT_NAMES::ATT_VERTNORMAL):
                    return vcg::tri::HasPerVertexNormal(_mesh);
            case(INT_ATT_NAMES::ATT_FACENORMAL):
                return vcg::tri::HasPerFaceNormal(_mesh);
            case(INT_ATT_NAMES::ATT_VERTCOLOR):
                return vcg::tri::HasPerVertexColor(_mesh);
            case(INT_ATT_NAMES::ATT_FACECOLOR):
                return vcg::tri::HasPerFaceColor(_mesh);
            /*case(INT_ATT_NAMES::ATT_FIXEDCOLOR):
                return true;*/
            case(INT_ATT_NAMES::ATT_VERTTEXTURE):
                return vcg::tri::HasPerVertexTexCoord(_mesh);
            case(INT_ATT_NAMES::ATT_WEDGETEXTURE):
                return vcg::tri::HasPerWedgeTexCoord(_mesh);
            case(INT_ATT_NAMES::ATT_VERTINDICES):
                return (_mesh.VN() != 0) && (_mesh.FN() != 0);
            case(INT_ATT_NAMES::ATT_EDGEINDICES):
                return vcg::tri::HasPerVertexFlags(_mesh) || ((_mesh.VN() != 0) && (_mesh.FN() == 0) && (_mesh.EN() == 0));
            default:
                return false;
            }
            return false;
        }

        bool checkBuffersAllocationStatus(InternalRendAtts& tobeallocated,InternalRendAtts& tobedeallocated,InternalRendAtts& tobeupdated) const
        {
            bool somethingtodo = false;
            tobedeallocated.reset();
            tobedeallocated.reset();
            tobeupdated.reset();
            
            //bool thereisreplicatedview = isThereAReplicatedPipelineView();
            InternalRendAtts meaningfulrequiredbyatleastoneview;
            InternalRendAtts probabilyuseless;
         
            for(typename ViewsMap::const_iterator it = _perviewreqatts.begin();it != _perviewreqatts.end();++it)
            {
                for(PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY;pm = next(pm))
                {
                    //If a primitive_modality is not rendered (== no att_VERTPOSITION) all the referred attributes by this view can be eventually deallocated IF they are not used
                    //by some other rendered primitive
                    //the vertindices is, as usual, a different case
                    if (it->second._intatts[size_t(pm)][INT_ATT_NAMES::ATT_VERTPOSITION])
                        meaningfulrequiredbyatleastoneview = InternalRendAtts::unionSet(meaningfulrequiredbyatleastoneview,it->second._intatts[size_t(pm)]);
                    else
                        probabilyuseless = InternalRendAtts::unionSet(probabilyuseless,it->second._intatts[size_t(pm)]);
                }
            }
            bool thereisreplicatedview = InternalRendAtts::replicatedPipelineNeeded(meaningfulrequiredbyatleastoneview);
            meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_VERTINDICES] &= !thereisreplicatedview;
            
            InternalRendAtts reallyuseless = InternalRendAtts::complementSet(probabilyuseless,meaningfulrequiredbyatleastoneview);
            
            bool switchreplicatedindexed = (!InternalRendAtts::replicatedPipelineNeeded(_currallocatedboatt) && thereisreplicatedview) || (InternalRendAtts::replicatedPipelineNeeded(_currallocatedboatt) && !thereisreplicatedview); 

            /*in some way the vertices number changed. If i use the indexed pipeline i have to deallocate/allocate/update the vertex indices*/
            bool numvertchanged = boExpectedSize(INT_ATT_NAMES::ATT_VERTPOSITION,thereisreplicatedview) != _bo[INT_ATT_NAMES::ATT_VERTPOSITION]->_size;
            bool vertindforcedupdate = numvertchanged && meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_VERTINDICES];

            InternalRendAtts probablytoallocate = InternalRendAtts::complementSet(meaningfulrequiredbyatleastoneview,_currallocatedboatt);
            InternalRendAtts probablytodeallocate = InternalRendAtts::complementSet(_currallocatedboatt,meaningfulrequiredbyatleastoneview);
            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                INT_ATT_NAMES boname(ii);
                if (_bo[boname] != NULL) 
                {
                    bool hasmeshattribute = hasMeshAttribute(boname);
                    bool isvalid = (_bo[boname]->_isvalid);               
                    bool notempty = (_bo[boname]->_size > 0);
                    if (boname != INT_ATT_NAMES::ATT_EDGEINDICES)
                    {
                        size_t sz = boExpectedSize(boname,thereisreplicatedview);

                        tobedeallocated[boname] =   (notempty && !hasmeshattribute) || 
                                                    (notempty && probablytodeallocate[boname]) ||
                                                    (notempty && reallyuseless[boname]) ||
                                                    (notempty && (_bo[boname]->_size != sz) && meaningfulrequiredbyatleastoneview[boname]) ||
                                                    (notempty && (boname == INT_ATT_NAMES::ATT_VERTINDICES) && (vertindforcedupdate));
                        tobeallocated[boname] = (hasmeshattribute && (sz > 0) && (sz != _bo[boname]->_size) && meaningfulrequiredbyatleastoneview[boname]) || 
                                                (hasmeshattribute && (sz > 0) && probablytoallocate[boname]) ||
                                                (hasmeshattribute && (boname == INT_ATT_NAMES::ATT_VERTINDICES) && (vertindforcedupdate));
                        tobeupdated[boname] = tobeallocated[boname] || (hasmeshattribute && (sz > 0) && !(isvalid) && meaningfulrequiredbyatleastoneview[boname]);
                    }
                    else 
                    {
                        bool meshchanged = ((_mesh.FN() != _meshtriangleswhenedgeindiceswerecomputed) || (_mesh.VN() != _meshverticeswhenedgeindiceswerecomputed));
                        tobedeallocated[INT_ATT_NAMES::ATT_EDGEINDICES] =   (notempty && !hasmeshattribute) || 
                                                                            (notempty && !meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_EDGEINDICES]) || 
                                                                            (notempty && !(isvalid) && meshchanged);
                        tobeallocated[INT_ATT_NAMES::ATT_EDGEINDICES] = (hasmeshattribute && meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_EDGEINDICES] && !(isvalid) && (meshchanged)) ||
                                                                        (hasmeshattribute && meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_EDGEINDICES] && !(isvalid) && !(_currallocatedboatt[INT_ATT_NAMES::ATT_EDGEINDICES]));
                        tobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES] = tobeallocated[INT_ATT_NAMES::ATT_EDGEINDICES] || 
                                                                      (hasmeshattribute && !(isvalid) && meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_EDGEINDICES]) ||
                                                                      (hasmeshattribute && switchreplicatedindexed && meaningfulrequiredbyatleastoneview[INT_ATT_NAMES::ATT_EDGEINDICES]);
                    } 
                }
                somethingtodo = somethingtodo || tobeallocated[boname] || tobedeallocated[boname] || tobeupdated[boname];
            }
            return !(somethingtodo);
        }

        bool manageAndFeedBuffersIfNeeded(const InternalRendAtts& tobeallocated,const InternalRendAtts& tobedeallocated,const InternalRendAtts& tobeupdated)
        {
            if (tobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES])
                updateEdgeVertIndVector();

            bool immediatemode = !(buffersMemoryManagementFunction(tobeallocated,tobedeallocated,tobeupdated));
            bool replicated = isThereAReplicatedPipelineView();

            if (immediatemode)
                return false;

            bool somethingtoupdate = false;
            for(unsigned int hh = 0;hh < INT_ATT_NAMES::enumArity();++hh)
                somethingtoupdate = somethingtoupdate || tobeupdated[hh];

            if (somethingtoupdate)
            {
                if (replicated)
                {
                    InternalRendAtts attributestobeupdated(tobeupdated);
                    //WARNING!In case we have to update the wedgetexture bo maybe (not always!) we must update also the other buffer already in memory
                    //cause the wedgetexture pipeline force a change in the order of the triangles in GPU.
                    //they are now ordered by the texture seam and not more by the triangle index!
                    if (tobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                        attributestobeupdated = _currallocatedboatt;
                    updateBuffersReplicatedPipeline(attributestobeupdated);
                }
                else
                    updateBuffersIndexedPipeline(tobeupdated);
                glFinish();
            }
            return true;
        }

        void updateEdgeVertIndVector()
        {
            _edge.clear();
            fillUniqueEdgeVector(_mesh,_edge);
            _meshverticeswhenedgeindiceswerecomputed = _mesh.VN();
            _meshtriangleswhenedgeindiceswerecomputed = _mesh.FN();
        }

        bool buffersMemoryManagementFunction(const InternalRendAtts& tobeallocated,const InternalRendAtts& tobedeallocated,const InternalRendAtts& tobeupdated)
        {
            //GLenum err = glGetError();
            bool replicated = isThereAReplicatedPipelineView();
            std::ptrdiff_t newallocatedmem = bufferObjectsMemoryRequired(tobeallocated);
            std::ptrdiff_t deallocatedmem = bufferObjectsMemoryRequired(tobedeallocated);
            ptrdiff_t zero = 0;
            std::ptrdiff_t changedsize = std::max(zero,newallocatedmem - deallocatedmem);
            //std::ptrdiff_t bomemoryrequiredbymesh = bufferObjectsMemoryRequired(_currallocatedboatt) - deallocatedmem + newallocatedmem;
            unsigned int ii = 0;
            for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
            {
                INT_ATT_NAMES boname(ii);
                //size_t sz = boExpectedSize(boname,replicated);
                //size_t dim = boExpectedDimension(boname,replicated);

                if (tobedeallocated[boname])
                    bufferDeAllocationRequested(boname); 
                ++ii;
            }

            if (!_gpumeminfo.isAdditionalMemoryAvailable(changedsize))
            {
                std::cout << "no additional memory available!!! memory required: " << changedsize << std::endl;
                ii = 0;
                for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
                {
                    INT_ATT_NAMES boname(ii);
                    size_t sz(boExpectedSize(boname,replicated));
                    //there are already valid mesh attributes properly allocated in memory but there is not enough gpu memory for the remaining mesh.
                    //we have to deallocate the previously allocated mesh attributes
                    if ((*it != NULL) && ((sz == (*it)->_size)))
                        bufferDeAllocationRequested(boname); 
                    ++ii;
                }
                _borendering = false;
                return false;
            }
            else
            {
                bool failedallocation = false;
                unsigned int ii = 0;
                typename std::vector<GLBufferObject*>::iterator it = _bo.begin();
                while((it != _bo.end()) && (!failedallocation))
                {
                    INT_ATT_NAMES boname(ii);
                    GLBufferObject* cbo = _bo[ii];
                    if (tobeallocated[boname])
                    {
                        cbo->_size = boExpectedSize(boname,replicated);
                        std::ptrdiff_t dim = boExpectedDimension(boname,replicated);
                        glGenBuffers(1, &cbo->_bohandle);
                        glBindBuffer(cbo->_target, cbo->_bohandle);
                        //we call glGetError BEFORE the glBufferData function in order to clean the error flag
                        GLenum err = glGetError();
                        //assert(err == GL_NO_ERROR);
                        glBufferData(cbo->_target, dim, NULL, GL_STATIC_DRAW);
                        err = glGetError();
                        //even if there according the MemoryInfo subclass there is enough space we were not able to allocate an attribute buffer object. We have to deallocate all the bos related to this mesh
                        failedallocation = (err == GL_OUT_OF_MEMORY) || (!_gpumeminfo.isAdditionalMemoryAvailable(dim));
                        if (!failedallocation)
                        {
                            //setBufferPointerEnableClientState(boname);
                            setBufferPointer(boname);
                            _gpumeminfo.acquiredMemory(dim);
                        }
                        cbo->_isvalid = !failedallocation;
                        _borendering = !failedallocation;
                        glBindBuffer(cbo->_target, 0);
                        _currallocatedboatt[boname] = !failedallocation;
                    }
                    else
                    {
                        //the arity of the attribute contained in the bo didn't change so i can use the old space without reallocating it
                        if (cbo != NULL)
                            cbo->_isvalid = cbo->_isvalid || tobeupdated[boname];
                    }
                    ++it;
                    ++ii;
                }
                if (failedallocation)
                   buffersDeAllocationRequested(_currallocatedboatt);
                _borendering = !failedallocation;
            }
            return _borendering;
        }

        bool updateBuffersIndexedPipeline(const InternalRendAtts& attributestobeupdated)
        {
            _chunkmap.clear();
            size_t vn = _mesh.VN();
            size_t tn = _mesh.FN();

            size_t facechunk = std::min(size_t(tn),_perbatchprim);
            size_t vertexchunk = std::min(size_t(vn),_perbatchprim);

            std::vector<vcg::Point3f> pv; //position vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                pv.resize(vertexchunk);

            std::vector<vcg::Point3f> nv; //per vertex normal vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                nv.resize(vertexchunk);

            std::vector<vcg::Color4b> cv; // Per vertex color vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR])
                cv.resize(vertexchunk);

            std::vector<float> tv; // per vertex texture coord vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                tv.resize(vertexchunk * 2);

            size_t chunkingpu = 0;

            for(size_t i=0;i<vn;++i)
            {
                size_t chunkindex = i % vertexchunk;
                if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                    pv[chunkindex].Import(_mesh.vert[i].cP());

                if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                {
                    nv[chunkindex].Import(_mesh.vert[i].cN());
                    nv[chunkindex].Normalize();
                }

                if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR])
                {
                    cv[chunkindex] = _mesh.vert[i].cC();
                }
                if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                {
                    tv[chunkindex*2+0] = _mesh.vert[i].cT().U();
                    tv[chunkindex*2+1] = _mesh.vert[i].cT().V();
                }

                if((i == vn - 1) || (chunkindex == vertexchunk - 1))
                {
                    size_t chunksize = vertexchunk;
                    if (i == vn - 1)
                        chunksize = chunkindex + 1;

                    if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                    {
                        GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTPOSITION];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&pv[0]);
                        //std::vector<vcg::Point3f> tmppv; //position vector
                        //if (attributestobeupdated[GLMeshAttributesInfo::ATT_VERTPOSITION])
                        //    tmppv.resize(vertexchunk);
                        //glGetBufferSubData(GL_ARRAY_BUFFER,0,buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&tmppv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);

                    }
                    if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                    {
                        GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTNORMAL];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&nv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR])
                    {
                        GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTCOLOR];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&cv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                    {
                        GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTTEXTURE];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&tv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    glFinish();
                    ++chunkingpu;
                }
            }

            pv.clear();
            nv.clear();
            cv.clear();
            tv.clear();

            chunkingpu = 0;
            std::vector<GLuint> ti(facechunk * 3);
            for(size_t i=0;i<tn;++i)
            {
                size_t chunkindex = i % facechunk;

                ti[chunkindex * 3 + 0] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(0)));
                ti[chunkindex * 3 + 1] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(1)));
                ti[chunkindex * 3 + 2] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(2)));

                if((i == tn - 1) || (chunkindex == facechunk - 1))
                {
                    size_t chunksize = facechunk;
                    if (i == tn - 1)
                        chunksize = chunkindex + 1;

                    if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTINDICES])
                    {
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[INT_ATT_NAMES::ATT_VERTINDICES]->_bohandle);
                        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,chunkingpu * facechunk *  _bo[INT_ATT_NAMES::ATT_VERTINDICES]->_components *  _bo[INT_ATT_NAMES::ATT_VERTINDICES]->getSizeOfGLType(),_bo[INT_ATT_NAMES::ATT_VERTINDICES]->_components *  _bo[INT_ATT_NAMES::ATT_VERTINDICES]->getSizeOfGLType() * chunksize,&ti[0]);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                    }
                    ++chunkingpu;
                }
            }
            if ((attributestobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES]) && (_edge.size() > 0))
            {
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[INT_ATT_NAMES::ATT_EDGEINDICES]->_bohandle);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,0,_bo[INT_ATT_NAMES::ATT_EDGEINDICES]->_components *  _edge.size() * _bo[INT_ATT_NAMES::ATT_EDGEINDICES]->getSizeOfGLType(),&_edge[0]);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            }
            return true;
        }

        bool isThereAReplicatedPipelineView() const
        {
            bool replicated = false;
            for(typename ViewsMap::const_iterator it = _perviewreqatts.begin();it != _perviewreqatts.end();++it)
            {
                //There is a replicated pipeline only if the att[ATT_VERTPOSITION] is true, otherwise is a spurious replicated view
                for(size_t pm = 0; pm < size_t(PR_ARITY);++pm) 
                    replicated = replicated || (InternalRendAtts::replicatedPipelineNeeded(it->second._intatts[pm]) && (it->second._pmmask.test(pm)));
            }
            return replicated;
        }

        bool isThereAnEdgesView() const
        {
            bool isthereaquadview = false;
            for(typename ViewsMap::const_iterator it = _perviewreqatts.begin();it != _perviewreqatts.end();++it)
                isthereaquadview = (it->second._intatts[size_t(PR_WIREFRAME_EDGES)][INT_ATT_NAMES::ATT_VERTPOSITION]) || isthereaquadview;
			return isthereaquadview;
        }


        bool updateBuffersReplicatedPipeline(const InternalRendAtts& attributestobeupdated)
        {
            size_t tn = _mesh.fn;

            size_t facechunk = std::min(size_t(tn),_perbatchprim);

            std::vector<vcg::Point3f> rpv; //position vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                rpv.resize(facechunk * 3);

            std::vector<vcg::Point3f> rnv; //per vertex normal vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                rnv.resize(facechunk * 3);

            std::vector<vcg::Point3f> rfnv; //per face normal vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_FACENORMAL])
                rfnv.resize(facechunk * 3);

            std::vector<vcg::Color4b> rcv; // Per vertex color vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR])
                rcv.resize(facechunk * 3);

            std::vector<vcg::Color4b> rfcv; // Per vertex color vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_FACECOLOR])
                rfcv.resize(facechunk * 3);

            std::vector<float> rtv; // per vertex texture coord vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                rtv.resize(facechunk * 3 * 2);

            std::vector<float> rwtv; // per wedge texture coord vector
            if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                rwtv.resize(facechunk * 3 * 2);

            size_t chunkingpu = 0;

            //it's a map containing for each texture seams n a vector of all the triangle index ranges having n has texture seam
            //Suppose that in a mesh we have
            //TXS_0{t0,t1,t2,t3}, TXS_4{t4,t5},TXS_0{t6},TXS_-1{t7,t8,t9},TXS_4{t10,t11}
            //so chunkMap will contain
            // -1 -> [<t7,t9>]
            //  0 -> [<t0,t3>,<t6,t6>]
            //  4 -> [<t4,t5>,<t10,t11>]
            //
            //if the map has no-texture coords at all in order to unify the code we fill the ChunkMap with texture seam -1 and a single triangle range going from face_0 to face_n-1


            if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE] || attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
            {
                _chunkmap.clear();
                if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                    fillchunkMap();
                else
                    if(attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                        _chunkmap[0].push_back(std::make_pair(0,tn-1));
            }

            //default case: no texture is required to be rendered but a non texture attribute has to be updated
            //we have to init the _chunkmap with just one entry (-1...that means no texture) referring all the triangles in the mesh
            if ((!_currallocatedboatt[INT_ATT_NAMES::ATT_VERTTEXTURE] && !_currallocatedboatt[INT_ATT_NAMES::ATT_WEDGETEXTURE]) &&
                (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION] ||
                attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL] || attributestobeupdated[INT_ATT_NAMES::ATT_FACENORMAL] ||
                attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR] || attributestobeupdated[INT_ATT_NAMES::ATT_FACECOLOR]))
            {
                _chunkmap.clear();
                _chunkmap[-1].push_back(std::make_pair(0,tn-1));
            }

            int t = 0;
            if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE] || attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
            {
                _texindnumtriangles.clear();
                _texindnumtriangles.resize(_chunkmap.size());
            }
            
            std::vector<GLuint> vpatlas;
            if (attributestobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES])
                vpatlas.resize(_mesh.VN(),UINT_MAX);

            int faceind = 0;
            size_t chunkindex = faceind;
			GLuint triangles = 0;


            for(ChunkMap::const_iterator mit = _chunkmap.begin();mit != _chunkmap.end();++mit)
            {
                for (ChunkVector::const_iterator cit = mit->second.begin();cit != mit->second.end();++cit)
                {
                    for(size_t indf = cit->first;indf<=cit->second;++indf)
                    {
                        chunkindex = faceind % facechunk;
                        if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                        {
                            rpv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->P());
                            rpv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->P());
                            rpv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->P());
                        }
                        if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                        {
                            rnv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->N().Normalize());
                            rnv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->N().Normalize());
                            rnv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->N().Normalize());
                        }

                        if (attributestobeupdated[INT_ATT_NAMES::ATT_FACENORMAL])
                        {
                            rfnv[chunkindex*3+0].Import(_mesh.face[indf].N().Normalize());
                            rfnv[chunkindex*3+1].Import(_mesh.face[indf].N().Normalize());
                            rfnv[chunkindex*3+2].Import(_mesh.face[indf].N().Normalize());
                        }

                        if ((attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR]))
                        {
                            rcv[chunkindex*3+0] = _mesh.face[indf].V(0)->C();
                            rcv[chunkindex*3+1] = _mesh.face[indf].V(1)->C();
                            rcv[chunkindex*3+2] = _mesh.face[indf].V(2)->C();
                        }

                        if ((attributestobeupdated[INT_ATT_NAMES::ATT_FACECOLOR]))
                        {
                            rfcv[chunkindex*3+0] = _mesh.face[indf].C();
                            rfcv[chunkindex*3+1] = _mesh.face[indf].C();
                            rfcv[chunkindex*3+2] = _mesh.face[indf].C();
                        }

                        if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                        {
                            rtv[chunkindex*6+0]=float(_mesh.face[indf].V(0)->T().U());
                            rtv[chunkindex*6+1]=float(_mesh.face[indf].V(0)->T().V());
                            rtv[chunkindex*6+2]=float(_mesh.face[indf].V(1)->T().U());
                            rtv[chunkindex*6+3]=float(_mesh.face[indf].V(1)->T().V());
                            rtv[chunkindex*6+4]=float(_mesh.face[indf].V(2)->T().U());
                            rtv[chunkindex*6+5]=float(_mesh.face[indf].V(2)->T().V());
                        }

                        if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                        {
                            rwtv[chunkindex*6+0]=float(_mesh.face[indf].WT(0).U());
                            rwtv[chunkindex*6+1]=float(_mesh.face[indf].WT(0).V());
                            rwtv[chunkindex*6+2]=float(_mesh.face[indf].WT(1).U());
                            rwtv[chunkindex*6+3]=float(_mesh.face[indf].WT(1).V());
                            rwtv[chunkindex*6+4]=float(_mesh.face[indf].WT(2).U());
                            rwtv[chunkindex*6+5]=float(_mesh.face[indf].WT(2).V());
                        }

                        if (attributestobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES])
                        {
                            for (int ii = 0;ii < 3;++ii)
                            {
                                size_t v = vcg::tri::Index(_mesh,_mesh.face[indf].V(ii));
                                if (vpatlas[v] == UINT_MAX)
                                    vpatlas[v] = faceind*3+ii;
                            }
                        }


                        if((faceind == tn - 1) || (chunkindex == facechunk - 1))
                        {
                            size_t chunksize = facechunk;
                            if (faceind == tn - 1)
                                chunksize = chunkindex + 1;

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTPOSITION])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTPOSITION];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rpv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTNORMAL])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTNORMAL];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rnv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_FACENORMAL])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_FACENORMAL];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rfnv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTCOLOR])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTCOLOR];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rcv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_FACECOLOR])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_FACECOLOR];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rfcv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_VERTTEXTURE];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rtv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                            {
                                GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_WEDGETEXTURE];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rwtv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            ++chunkingpu;
                        }
                        ++faceind;
                    }
                    triangles += cit->second - cit->first + 1;
                }

                if ((attributestobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES]) && (_edge.size() > 0)) 
                {
                    for(typename std::vector<EdgeVertInd>::iterator it = _edge.begin();it != _edge.end();++it)
                    {
                        it->_v[0] = vpatlas[it->_v[0]];
                        it->_v[1] = vpatlas[it->_v[1]];
                    }

                    GLBufferObject* buffobj = _bo[INT_ATT_NAMES::ATT_EDGEINDICES];
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffobj->_bohandle);
                    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,0,buffobj->_components * buffobj->getSizeOfGLType() * _edge.size(),&_edge[0]);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                }

                if (attributestobeupdated[INT_ATT_NAMES::ATT_WEDGETEXTURE] || attributestobeupdated[INT_ATT_NAMES::ATT_VERTTEXTURE])
                    _texindnumtriangles[t] = std::make_pair(mit->first,triangles);
                ++t;
            }

            //return (k != tn)
            //    throw MeshLabException("Mesh has not been properly partitioned");
            return true;
        }

        void buffersDeAllocationRequested(const InternalRendAtts& rq)
        {
            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                INT_ATT_NAMES boname(ii);
                if ((_bo[ii] != NULL) && (rq[boname]))
                    bufferDeAllocationRequested(boname);
            }
        }

        void bufferDeAllocationRequested(INT_ATT_NAMES att)
        {
            unsigned int ind(att);
            if ((ind < 0) || (ind >= (unsigned int) _bo.size()))
                return;
            GLBufferObject* bobj = _bo[ind];
            if (bobj == NULL)
                return;

            if ((att != INT_ATT_NAMES::ATT_VERTINDICES) && (att != INT_ATT_NAMES::ATT_EDGEINDICES) /*&& (att != INT_ATT_NAMES::ATT_FIXEDCOLOR)*/)
            {
                glDisableClientState(bobj->_clientstatetag);
            }
			//glBufferData(bobj->_target, sizeof(vcg::Point3f)*_primitivebatch, 0, GL_DYNAMIC_DRAW);
            glDeleteBuffers(1,&(bobj->_bohandle));
            glFlush();
            glFinish();

            if (bobj->_size > 0)
                //we don't use dim cause dim is the value that is going to be allocated, instead use (*it)->_size * (*it)->getSizeOfGLType() is the value already in the buffer
                _gpumeminfo.releasedMemory(bobj->_size * bobj->getSizeOfGLType());
            bobj->_isvalid = false;
            bobj->_size = 0;
            _currallocatedboatt[att] = false; 
        }

        std::ptrdiff_t bufferObjectsMemoryRequired(const InternalRendAtts& rqatt) const
        {
            bool replicated = InternalRendAtts::replicatedPipelineNeeded(rqatt);
            std::ptrdiff_t result(0);

            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                INT_ATT_NAMES nm(ii);
                if (rqatt[nm])
                    result += (std::ptrdiff_t) boExpectedDimension(nm,replicated);
            }
            return result;
        }

        //expected number of cells should have the required bo
        //generateindex is true when i have a triangle based mesh
        //				is false when i have a point based mesh
        size_t boExpectedSize(INT_ATT_NAMES name,bool replicatedpipeline) const
        {
            try
            {
                
                GLBufferObject& cbo = *_bo.at((unsigned int) name);
                size_t vertnum(_mesh.VN());
                size_t facenum(_mesh.FN());

                switch((unsigned int) name)
                {
                case(INT_ATT_NAMES::ATT_VERTPOSITION):
                case(INT_ATT_NAMES::ATT_VERTNORMAL):
                case(INT_ATT_NAMES::ATT_VERTCOLOR):
                case(INT_ATT_NAMES::ATT_VERTTEXTURE):
                    {
                        if (replicatedpipeline)
                            return facenum * 3 * cbo._components;
                        else
                            return vertnum * cbo._components;
                    }

                case(INT_ATT_NAMES::ATT_FACENORMAL):
                case(INT_ATT_NAMES::ATT_FACECOLOR):
                case(INT_ATT_NAMES::ATT_WEDGETEXTURE):
                    {
                        if (replicatedpipeline)
                            return facenum * 3 * cbo._components;
                        else
                            return 0;
                    }
                case(INT_ATT_NAMES::ATT_VERTINDICES):
                    {
                        if (replicatedpipeline)
                            return 0;
                        else
                            return facenum * cbo._components;
                    }
                case(INT_ATT_NAMES::ATT_EDGEINDICES):
                    {
                        return _edge.size() * cbo._components;
                    }

                default : break;
                }
            }
            catch(std::out_of_range& /*exc*/)
            {
                return 0;
            }
            return 0;
        }

        size_t boExpectedDimension(INT_ATT_NAMES name,bool replicatedpipeline) const
        {
            try
            {
                size_t sz = boExpectedSize(name,replicatedpipeline);
                unsigned int ind = (unsigned int) name;
                GLBufferObject* cbo = _bo.at(ind);
                if (cbo == NULL)
                    return 0;
                else
                    return sz * cbo->getSizeOfGLType();
            }
            catch(std::out_of_range& /*exc*/)
            {
                return 0;
            }
            return 0;
        }

		void drawFun(const PVData& dt, const std::vector<GLuint>& textid = std::vector<GLuint>()) const
		{
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glMultMatrix(_tr);

			if ((dt._glopts != NULL) && (dt._glopts->_perbbox_enabled))
				drawBBox(dt._glopts);

			if (dt.isPrimitiveActive(PR_SOLID))
			{
				bool somethingmore = dt.isPrimitiveActive(PR_WIREFRAME_EDGES) || dt.isPrimitiveActive(PR_WIREFRAME_TRIANGLES) || dt.isPrimitiveActive(PR_POINTS);
				if (somethingmore)
				{
					glEnable(GL_POLYGON_OFFSET_FILL);
					glPolygonOffset(1.0, 1);
				}
				drawFilledTriangles(dt._intatts[size_t(PR_SOLID)], dt._glopts, textid);
				if (somethingmore)
					glDisable(GL_POLYGON_OFFSET_FILL);
			}

			if (dt.isPrimitiveActive(PR_WIREFRAME_EDGES) || dt.isPrimitiveActive(PR_WIREFRAME_TRIANGLES))
			{
				//InternalRendAtts tmpatts = atts;
				bool pointstoo = dt.isPrimitiveActive(PR_POINTS);

				if (pointstoo)
				{
					glEnable(GL_POLYGON_OFFSET_FILL);
					glPolygonOffset(1.0, 1);
				}
				bool solidtoo = dt.isPrimitiveActive(PR_SOLID);

				/*EDGE    |     TRI     |    DRAW
				---------------------------------
				  TRUE         TRUE         EDGE
				  TRUE         FALSE        EDGE
				  FALSE        TRUE         TRI
				  FALSE        FALSE        NOTHING */

				if (dt.isPrimitiveActive(PR_WIREFRAME_EDGES))
					drawEdges(dt._intatts[size_t(PR_WIREFRAME_EDGES)], dt._glopts);
				else
				{
					if (dt.isPrimitiveActive(PR_WIREFRAME_TRIANGLES))
					{
						drawWiredTriangles(dt._intatts[size_t(PR_WIREFRAME_TRIANGLES)], dt._glopts, textid);
					}
				}

				if (pointstoo || solidtoo)
					glDisable(GL_POLYGON_OFFSET_FILL);
			}
			if (dt.isPrimitiveActive(PR_POINTS))
				drawPoints(dt._intatts[size_t(PR_POINTS)], dt._glopts,textid);

			glPopMatrix();
			glPopAttrib();
			glFlush();
			glFinish();
		}

        void drawFilledTriangles(const InternalRendAtts& req,const GL_OPTIONS_DERIVED_TYPE* glopts,const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const
        {
            if (_mesh.VN() == 0)
                return;

            glPushAttrib(GL_ALL_ATTRIB_BITS);
  
            bool isgloptsvalid = (glopts != NULL);
            
            if (isgloptsvalid && glopts->_persolid_noshading)
                glDisable(GL_LIGHTING);
            else
                if ((!isgloptsvalid) ||  (req[INT_ATT_NAMES::ATT_VERTNORMAL]) || (req[INT_ATT_NAMES::ATT_FACENORMAL]))
                {
                    glEnable(GL_LIGHTING);
                }

            glEnable(GL_COLOR_MATERIAL);
            if ((isgloptsvalid) && (glopts->_persolid_fixed_color_enabled))
                glColor(glopts->_persolid_fixed_color);
            else
            {
                if ((isgloptsvalid) && (glopts->_persolid_mesh_color_enabled))
                    glColor(_mesh.C());
                else
                {
                    if ((req[INT_ATT_NAMES::ATT_VERTCOLOR]) || (req[INT_ATT_NAMES::ATT_FACECOLOR]))
                        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
                    else
                        glColor(vcg::Color4b(vcg::Color4b::LightGray));
                }
            }

            if (isBORenderingAvailable())
                drawTrianglesBO(req,textureindex);
            else
                drawTrianglesIM(req,textureindex);

            glPopAttrib();
        }


        void drawWiredTriangles(const InternalRendAtts& req,const GL_OPTIONS_DERIVED_TYPE* glopts,const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const
        {
            if (_mesh.VN() == 0)
                return;
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            
           
            bool isgloptsvalid = (glopts != NULL);

            if (isgloptsvalid && glopts->_perwire_noshading)
                glDisable(GL_LIGHTING);
            else 
                if ((!isgloptsvalid) || (req[INT_ATT_NAMES::ATT_VERTNORMAL]))
                {
                    glEnable(GL_LIGHTING);
                }
                
            glEnable(GL_COLOR_MATERIAL);
            if ((isgloptsvalid) && (glopts->_perwire_fixed_color_enabled))
                glColor(glopts->_perwire_fixed_color);
            else
            {
                if ((isgloptsvalid) && (glopts->_perwire_mesh_color_enabled))
                    glColor(_mesh.C());
                else
                {
                    if (req[INT_ATT_NAMES::ATT_VERTCOLOR])
                        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
                    else
                        glColor(vcg::Color4b(vcg::Color4b::DarkGray));
                }
            }
            float linewidth = 1.0f;
            if (isgloptsvalid)
                linewidth = glopts->_perwire_wirewidth;
            glLineWidth(linewidth);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

            if (isBORenderingAvailable())
                drawTrianglesBO(req,textureindex);
            else
                drawTrianglesIM(req,textureindex);

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glPopAttrib();
        }

        void drawTrianglesBO(const InternalRendAtts& req,const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const
        {
            updateClientState(req);

            bool replicated = InternalRendAtts::replicatedPipelineNeeded(_currallocatedboatt);

            if (replicated)
            {
                //qDebug("Replicated drawing");
                int firsttriangleoffset = 0;
                if(!req[INT_ATT_NAMES::ATT_VERTTEXTURE] && !req[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                {
                    glDisable(GL_TEXTURE_2D);
                    glDrawArrays(GL_TRIANGLES,0,_mesh.fn * 3);
                }
                else
                {
                    glEnable(GL_TEXTURE_2D);
                    for(std::vector< std::pair<short,GLuint> >::const_iterator it = _texindnumtriangles.begin();it != _texindnumtriangles.end();++it)
                    {
                        if ((it->first != -1) && (it->first < textureindex.size()))
                            glBindTexture(GL_TEXTURE_2D,textureindex[it->first]);
                        else
                            glBindTexture(GL_TEXTURE_2D,0);
                        glDrawArrays(GL_TRIANGLES,firsttriangleoffset,it->second * 3 - firsttriangleoffset);
                        firsttriangleoffset = it->second * 3;
                    }
                    glBindTexture(GL_TEXTURE_2D,0);
                    glDisable(GL_TEXTURE_2D);
                }

            }
            else
            {
                if(req[INT_ATT_NAMES::ATT_VERTTEXTURE])
                {
                    if (textureindex.size() > 0)
                    {
                        glEnable(GL_TEXTURE_2D);
                        glBindTexture(GL_TEXTURE_2D,textureindex[0]);
                    }
                }
                else
                    glDisable(GL_TEXTURE_2D);

                if  (_bo[INT_ATT_NAMES::ATT_VERTINDICES]->_isvalid)
                {
                    //qDebug("Indexed drawing");
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,_bo[INT_ATT_NAMES::ATT_VERTINDICES]->_bohandle);
                    glDrawElements( GL_TRIANGLES, _mesh.FN() * _bo[INT_ATT_NAMES::ATT_VERTINDICES]->_components,GL_UNSIGNED_INT ,NULL);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);   

                }

                glBindTexture(GL_TEXTURE_2D,0);
                glDisable(GL_TEXTURE_2D);
            }
            InternalRendAtts tmp;
            updateClientState(tmp);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
            glBindBuffer(GL_ARRAY_BUFFER,0);
        }

        void drawTrianglesIM(const InternalRendAtts& req,const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const
        {
            if(_mesh.fn==0)
                return;

            bool vn = req[INT_ATT_NAMES::ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(_mesh);
            bool fn = req[INT_ATT_NAMES::ATT_FACENORMAL] && vcg::tri::HasPerFaceNormal(_mesh);
            bool vc = req[INT_ATT_NAMES::ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(_mesh);
            bool fc = req[INT_ATT_NAMES::ATT_FACECOLOR] && vcg::tri::HasPerFaceColor(_mesh);
            bool vt = req[INT_ATT_NAMES::ATT_VERTTEXTURE] && vcg::tri::HasPerVertexTexCoord(_mesh);
            bool wt = req[INT_ATT_NAMES::ATT_WEDGETEXTURE] && vcg::tri::HasPerWedgeTexCoord(_mesh);


            //typename MESHTYPE::FaceContainer::iterator fp;
            typename MESH_TYPE::FaceIterator fi = _mesh.face.begin();

            short curtexname=-1;
            if(wt)
            {
                curtexname=(*fi).WT(0).n();
                if ((curtexname >= 0) && (curtexname < (int)textureindex.size()))
                {
                    glEnable(GL_TEXTURE_2D);
                    glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
                }
                else
                {
                    glDisable(GL_TEXTURE_2D);
                }
            }

            if(vt && !textureindex.empty()) // in the case of per vertex tex coord we assume that we have a SINGLE texture.
            {
                curtexname = 0;
                glEnable(GL_TEXTURE_2D);
                glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
            }

            glBegin(GL_TRIANGLES);

            while(fi!=_mesh.face.end())
            {
                typename MESH_TYPE::FaceType & f = *fi;
                if(!f.IsD())
                {
                    if(wt)
                        if(f.WT(0).n() != curtexname)
                        {
                            curtexname=(*fi).WT(0).n();
                            glEnd();

                            if (curtexname >= 0)
                            {
                                glEnable(GL_TEXTURE_2D);
                                if(!textureindex.empty())
                                    glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
                            }
                            else
                            {
                                glDisable(GL_TEXTURE_2D);
                            }

                            glBegin(GL_TRIANGLES);
                        }

                        if(fn)
                            glNormal(f.cN());
                        if(vn)
                            glNormal(f.V(0)->cN());

                        if(fc)
                            glColor(f.C());
                        if(vc)
                            glColor(f.V(0)->C());
                        if(vt)
                            glTexCoord(f.V(0)->T().P());
                        if(wt)
                            glTexCoord(f.WT(0).t(0));
                        glVertex(f.V(0)->P());

                        if(vn)
                            glNormal(f.V(1)->cN());
                        if(vc)
                            glColor(f.V(1)->C());
                        if(vt)
                            glTexCoord(f.V(1)->T().P());
                        if(wt)
                            glTexCoord(f.WT(1).t(0));
                        glVertex(f.V(1)->P());

                        if(vn)
                            glNormal(f.V(2)->cN());
                        if(vc)
                            glColor(f.V(2)->C());
                        if(vt)
                            glTexCoord(f.V(2)->T().P());
                        if(wt)
                            glTexCoord(f.WT(2).t(0));
                        glVertex(f.V(2)->P());
                }
                ++fi;
            }

            glEnd();
        }

        void drawPoints(const InternalRendAtts& req,GL_OPTIONS_DERIVED_TYPE* glopts, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const
        {
            if (_mesh.VN() == 0)
                return;
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            

            bool isgloptsvalid = (glopts != NULL);        

            
            if ((isgloptsvalid && glopts->_perpoint_noshading) || (isgloptsvalid && glopts->_perpoint_dot_enabled))
                glDisable(GL_LIGHTING);
            else 
                if ((!isgloptsvalid) ||  req[INT_ATT_NAMES::ATT_VERTNORMAL])
                {
                    glEnable(GL_LIGHTING);
                }

            glEnable(GL_COLOR_MATERIAL);
            if ((isgloptsvalid) && ((glopts->_perpoint_fixed_color_enabled) || (glopts->_perpoint_mesh_color_enabled))){
                if (glopts->_perpoint_fixed_color_enabled)
                    glColor(glopts->_perpoint_fixed_color);
                else
                    glColor(_mesh.C());
            }

            if (req[INT_ATT_NAMES::ATT_VERTCOLOR])       
                glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
         

			if (req[INT_ATT_NAMES::ATT_VERTTEXTURE])
			{
				glEnable(GL_TEXTURE_2D);
				if (textureindex.size() > 0)
					glBindTexture(GL_TEXTURE_2D, textureindex[0]);
				else
					glBindTexture(GL_TEXTURE_2D, 0);
			}
			else
				glDisable(GL_TEXTURE_2D);
            //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[GLMeshAttributesInfo::ATT_VERTINDEX]->_bohandle);
           
            if (glopts != NULL)
            {
				if (!glopts->_perpoint_dot_enabled)
					glPointSize(glopts->_perpoint_pointsize);
                if ((glopts->_perpoint_pointsmooth_enabled) || (glopts->_perpoint_dot_enabled))
                    glEnable(GL_POINT_SMOOTH);
                else 
                    glDisable(GL_POINT_SMOOTH);
                if(glopts->_perpoint_pointattenuation_enabled)
                {
                    vcg::Matrix44<typename MESH_TYPE::ScalarType> mat;
                    glGetv(GL_MODELVIEW_MATRIX,mat);
                    vcg::Point3<typename MESH_TYPE::ScalarType> c =_mesh.bbox.Center();
                    float camDist = (float)Norm(mat*c);
                    float quadratic[] = { 0.0f, 0.0f, 1.0f/(camDist*camDist) , 0.0f };
                    glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
                    glPointParameterf( GL_POINT_SIZE_MAX, 16.0f );
                    glPointParameterf( GL_POINT_SIZE_MIN, 1.0f );
                }
                else
                {
                    float quadratic[] = { 1.0f, 0.0f, 0.0f};
                    glPointParameterfv( GL_POINT_DISTANCE_ATTENUATION, quadratic );
                    float pointsize = 1.0f;
                    if (isgloptsvalid)
                        pointsize = glopts->_perpoint_pointsize;
                    glPointSize(pointsize);
                }
				if (glopts->_perpoint_dot_enabled)
				{
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glDepthRange(0.0, 0.9999);
					glDepthFunc(GL_LEQUAL);
					glPointSize(glopts->_perpoint_pointsize + 0.5);
				}
            }
            if (isBORenderingAvailable())
                drawPointsBO(req);
            else
                drawPointsIM(req);

			if ((glopts != NULL) && (glopts->_perpoint_dot_enabled))
			{
				float psize = 0.0001f;
				if ((glopts->_perpoint_pointsize - 1) > 0)
					psize = (glopts->_perpoint_pointsize - 1);
				glPointSize(psize);
				if (isBORenderingAvailable())
					drawPointsBO(req);
				else
					drawPointsIM(req);
			}
            glPopAttrib();
        }

        void drawPointsBO(const InternalRendAtts& req) const
        {
            size_t pointsnum = _mesh.VN();
            if (InternalRendAtts::replicatedPipelineNeeded(_currallocatedboatt))
                pointsnum = _mesh.FN() * 3;
            updateClientState(req);
            glDrawArrays(GL_POINTS,0,GLsizei(pointsnum));
            /*disable all client state buffers*/
            InternalRendAtts tmp;
            updateClientState(tmp);
        }

        void drawPointsIM(const InternalRendAtts& req) const
        {
            bool vn = req[INT_ATT_NAMES::ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(_mesh);
            bool vc = req[INT_ATT_NAMES::ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(_mesh);
            bool vt = req[INT_ATT_NAMES::ATT_VERTTEXTURE] && vcg::tri::HasPerVertexTexCoord(_mesh);
            

            glBegin(GL_POINTS);
            for(typename MESH_TYPE::VertexIterator vi=_mesh.vert.begin();vi!=_mesh.vert.end();++vi)
            {
                if(!(*vi).IsD())
                {
                    if(vn) glNormal((*vi).cN());
                    if(vc) glColor((*vi).C());
                    if(vt) glTexCoord((*vi).T().P());
                    glVertex((*vi).P());
                }
            }
            glEnd();
        }

        void drawEdges(const InternalRendAtts& req,GL_OPTIONS_DERIVED_TYPE* glopts) const
        {
            if (_mesh.VN() == 0)
                return;
            glPushAttrib(GL_ALL_ATTRIB_BITS);
           
            bool isgloptsvalid = (glopts != NULL);

            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

            if (isgloptsvalid && glopts->_perwire_noshading)
                glDisable(GL_LIGHTING);
            else 
                if ((!isgloptsvalid) || (req[INT_ATT_NAMES::ATT_VERTNORMAL]))
                {
                    glEnable(GL_LIGHTING);
                }

            bool colordefinedenabled = (isgloptsvalid) && ((glopts->_perwire_fixed_color_enabled) || (glopts->_perwire_mesh_color_enabled));

            if (!(isgloptsvalid) || colordefinedenabled)
            {
                vcg::Color4b tmpcol = vcg::Color4b(vcg::Color4b::DarkGray);
				if (colordefinedenabled)
				{
					if (glopts->_perwire_fixed_color_enabled)
						tmpcol = glopts->_perwire_fixed_color;
					else
						tmpcol = _mesh.cC();
				}
                glColor(tmpcol);
            }

            glDisable(GL_TEXTURE_2D);
            //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[GLMeshAttributesInfo::ATT_VERTINDEX]->_bohandle);

            float linewidth = 1.0f;
            if (isgloptsvalid)
                linewidth = glopts->_perwire_wirewidth;
            glLineWidth(linewidth);

            if (isBORenderingAvailable())
                drawEdgesBO(req);
            else
                drawEdgesIM(req);
            //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

            /*disable all client state buffers*/
            glPopAttrib();
        }
        
        void drawEdgesBO(const InternalRendAtts& req) const
        {
            if  (_bo[INT_ATT_NAMES::ATT_EDGEINDICES]->_isvalid)
            {
                //qDebug("Indexed drawing");
                updateClientState(req);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,_bo[INT_ATT_NAMES::ATT_EDGEINDICES]->_bohandle);
                glDrawElements( GL_LINES, _edge.size() * _bo[INT_ATT_NAMES::ATT_EDGEINDICES]->_components,GL_UNSIGNED_INT ,NULL);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                InternalRendAtts tmp;
                updateClientState(tmp);
            }
        }

        void drawEdgesIM(const InternalRendAtts& req) const
        {
            typename MESH_TYPE::FaceIterator fi = _mesh.face.begin();

            bool vn = req[INT_ATT_NAMES::ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(_mesh);
            bool vc = req[INT_ATT_NAMES::ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(_mesh);

            glBegin(GL_LINES);

            while(fi!=_mesh.face.end())
            {
                typename MESH_TYPE::FaceType & f = *fi;

                if(!f.IsD())
                {
                    if (!f.IsF(0)) 
                    {
                        if(vn)	glNormal(f.V(0)->cN());
                        if(vc)	glColor(f.V(0)->C());
                        glVertex(f.V(0)->P());

                        if(vn)	glNormal(f.V(1)->cN());
               
                        if(vc)	glColor(f.V(1)->C());
                        glVertex(f.V(1)->P());
                    }

                    if (!f.IsF(1)) 
                    {
                        if(vn)	glNormal(f.V(1)->cN());
                        if(vc)	glColor(f.V(1)->C());
                        glVertex(f.V(1)->P());

                        if(vn)	glNormal(f.V(2)->cN());
                        if(vc)	glColor(f.V(2)->C());
                        glVertex(f.V(2)->P());
                    }

                    if (!f.IsF(2)) 
                    {
                        if(vn)	glNormal(f.V(2)->cN());
                        if(vc)	glColor(f.V(2)->C());
                        glVertex(f.V(2)->P());

                        if(vn)	glNormal(f.V(0)->cN());
                        if(vc)	glColor(f.V(0)->C());
                        glVertex(f.V(0)->P());
                    }

                }
                ++fi;
            }

            glEnd();
        }

        void drawBBox(GL_OPTIONS_DERIVED_TYPE* glopts) const
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            bool isgloptsvalid = (glopts != NULL);

            glDisable(GL_LIGHTING);
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

            if ((isgloptsvalid) && (glopts->_perbbox_fixed_color_enabled))
                glColor(glopts->_perbbox_fixed_color);
            else
            {
                if ((isgloptsvalid) && (glopts->_perbbox_mesh_color_enabled))
                    glColor(_mesh.C());
                else
                    glColor(vcg::Color4b(vcg::Color4b::White));
            }
            if (isBORenderingAvailable())
                drawBBoxBO();
            else
                drawBBoxIM();
            glPopAttrib();
        }

        void drawBBoxBO() const
        {
            vcg::Box3<typename MESH_TYPE::ScalarType>& b = _mesh.bbox;

            GLuint bbhandle;
            glGenBuffers(1,&bbhandle);
            std::vector<vcg::Point3f> bbox(12 * 2);

            //0
            bbox[0] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
            bbox[1] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);

            //1
            bbox[2] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);
            bbox[3] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);

            //2
            bbox[4] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);
            bbox[5] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);

            //3
            bbox[6] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);
            bbox[7] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);

            //4
            bbox[8] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
            bbox[9] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);

            //5
            bbox[10] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);
            bbox[11] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);

            //6
            bbox[12] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);
            bbox[13] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);

            //7
            bbox[14] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);
            bbox[15] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);

            //8
            bbox[16] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);
            bbox[17] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);

            //9
            bbox[18] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);
            bbox[19] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);

            //10
            bbox[20] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
            bbox[21] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);

            //11
            bbox[22] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
            bbox[23] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);

            glBindBuffer(GL_ARRAY_BUFFER,bbhandle);
            glBufferData(GL_ARRAY_BUFFER, 12 * 2 * sizeof(vcg::Point3f), &(bbox[0]), GL_STATIC_DRAW);
            glVertexPointer(3,GL_FLOAT,0,0);
            glBindBuffer(GL_ARRAY_BUFFER,0);
            glEnableClientState(GL_VERTEX_ARRAY);
            glDrawArrays(GL_LINES,0,24);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDeleteBuffers(1,&bbhandle);
        }

        void drawBBoxIM() const
        {
            vcg::Box3<typename MESH_TYPE::ScalarType>& b = _mesh.bbox;

            glBegin(GL_LINE_STRIP);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
            glVertex3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);
            glVertex3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);
            glVertex3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
            glEnd();
            glBegin(GL_LINE_STRIP);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);
            glVertex3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);
            glVertex3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);
            glVertex3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);
            glEnd();
            glBegin(GL_LINES);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
            glVertex3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);

            glVertex3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);
            glVertex3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);

            glVertex3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);
            glVertex3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);

            glVertex3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);
            glVertex3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
            glEnd();
        }

        void updateClientState(const InternalRendAtts& req) const 
        {
            int ii = 0;
            for(typename std::vector<GLBufferObject*>::const_iterator it = _bo.begin();it != _bo.end();++it)
            {
                INT_ATT_NAMES boname(ii);
                if ((boname != INT_ATT_NAMES::ATT_VERTINDICES) && (boname != INT_ATT_NAMES::ATT_EDGEINDICES) /*&& (boname != INT_ATT_NAMES::ATT_FIXEDCOLOR)*/)
                {
                    if (req[boname] && _currallocatedboatt[boname] && (*it != NULL))
                    {
                        glBindBuffer((*it)->_target, (*it)->_bohandle);
                        setBufferPointer(boname);
                        glEnableClientState((*it)->_clientstatetag);
                        glBindBuffer((*it)->_target, 0);
                    }
                    else
                    {
                        glBindBuffer((*it)->_target, (*it)->_bohandle);
                        disableClientState(boname,req);
                        glBindBuffer((*it)->_target, 0);
                    }
                }
                ++ii;
            }
        }

        void setBufferPointer( INT_ATT_NAMES boname) const
        {
            unsigned int ii = boname;
            if (ii >= INT_ATT_NAMES::enumArity())
                return;
            GLBufferObject* cbo = _bo[ii];
            if (cbo == NULL)
                return;

            switch(ii)
            {
            case(INT_ATT_NAMES::ATT_VERTPOSITION):
                {
                    glVertexPointer(GLint(cbo->_components), cbo->_gltype, GLsizei(0), 0);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTNORMAL):
            case(INT_ATT_NAMES::ATT_FACENORMAL):
                {
                    glNormalPointer(cbo->_gltype, GLsizei(0), 0);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTCOLOR):
            case(INT_ATT_NAMES::ATT_FACECOLOR):
                {
                    glColorPointer(GLint(cbo->_components), cbo->_gltype, GLsizei(0), 0);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTTEXTURE):
            case(INT_ATT_NAMES::ATT_WEDGETEXTURE):
                {
                    glTexCoordPointer(GLint(cbo->_components), cbo->_gltype,GLsizei(0), 0);
                    break;
                }
            default : break;
            }
        }

        void disableClientState( INT_ATT_NAMES boname,const RendAtts& req) const
        {

            if (boname >= INT_ATT_NAMES::enumArity())
                return;

            switch(boname)
            {
            case(INT_ATT_NAMES::ATT_VERTPOSITION):
                {
                    glDisableClientState(GL_VERTEX_ARRAY);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTNORMAL):
            case(INT_ATT_NAMES::ATT_FACENORMAL):
                {
                    if (!req[INT_ATT_NAMES::ATT_VERTNORMAL] && !req[INT_ATT_NAMES::ATT_FACENORMAL])
                        glDisableClientState(GL_NORMAL_ARRAY);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTCOLOR):
            case(INT_ATT_NAMES::ATT_FACECOLOR):
                {
                    if (!req[INT_ATT_NAMES::ATT_VERTCOLOR] && !req[INT_ATT_NAMES::ATT_FACECOLOR])
                        glDisableClientState(GL_COLOR_ARRAY);
                    break;
                }
            case(INT_ATT_NAMES::ATT_VERTTEXTURE):
            case(INT_ATT_NAMES::ATT_WEDGETEXTURE):
                {
                    if (!req[INT_ATT_NAMES::ATT_VERTTEXTURE] && !req[INT_ATT_NAMES::ATT_WEDGETEXTURE])
                        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
                    break;
                }
            default:
                {
                    break;
                }

            }
        }

        void fillchunkMap()
        {
            if (!vcg::tri::HasPerWedgeTexCoord(_mesh))
                return;
            _chunkmap.clear();
            typename MESH_TYPE::FaceIterator infrange = _mesh.face.begin();
            short texind = std::numeric_limits<short>::max();
            int hh = 0;
            for(typename MESH_TYPE::FaceIterator fit = _mesh.face.begin();fit != _mesh.face.end();++fit)
            {
                if (fit->WT(0).N() != texind)
                {
                    if ((texind != std::numeric_limits<short>::max()) || (fit == _mesh.face.end() - 1))
                    {
                        GLuint lowind = std::distance(_mesh.face.begin(),infrange);
                        GLuint topind = std::distance(_mesh.face.begin(),fit) - 1;
                        _chunkmap[texind].push_back(std::make_pair(lowind,topind));
                        infrange = fit;
                    }
                    texind = fit->WT(0).N();
                }
                ++hh;
            }
            _chunkmap[texind].push_back(std::make_pair(std::distance(_mesh.face.begin(),infrange),std::distance(_mesh.face.begin(),_mesh.face.end() - 1)));
        }

        void debug(const InternalRendAtts& tobeallocated,const InternalRendAtts& tobedeallocated,const InternalRendAtts& tobeupdated)
        {
            _loginfo.reset();
            _loginfo._tobedeallocated = std::string("to_be_deallocated: ");
            _loginfo._tobeallocated = std::string("to_be_allocated: ");
            _loginfo._tobeupdated = std::string("to_be_updated: ");

            std::string truestring("true");
            std::string falsestring("false");
            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                std::string deallocres(falsestring);
                if (tobedeallocated[ii])
                    deallocres = truestring;
                _loginfo._tobedeallocated +=  deallocres + " ";

                std::string allocres(falsestring);
                if (tobeallocated[ii])
                    allocres = truestring;
                _loginfo._tobeallocated +=  allocres + " ";

                std::string upres(falsestring);
                if (tobeupdated[ii])
                    upres = truestring;
                _loginfo._tobeupdated +=  upres + " ";
            }

            _loginfo._tobedeallocated = std::string("[") + _loginfo._tobedeallocated + std::string("]");
            _loginfo._tobeallocated = std::string("[") + _loginfo._tobeallocated + std::string("]");
            _loginfo._tobeupdated = std::string("[") + _loginfo._tobeupdated + std::string("]");

         
            int hh = 0;
            
            for(typename ViewsMap::const_iterator it = _perviewreqatts.begin();it != _perviewreqatts.end();++it)
            {
                std::stringstream tmpstream;
                tmpstream << "view_" << hh << ":\n";
                for(size_t pm = 0; pm < size_t(PR_ARITY); ++pm)
                {
                    tmpstream << DebugInfo::primitiveName(pm) << " ";

                    for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
                    {
                        std::string res = falsestring;
                        if (it->second._intatts[pm][ii])
                            res = truestring;
                        tmpstream << "att[" << ii << "]=" << res << " ";
                    }
                    tmpstream << std::endl;
                }
                _loginfo._perviewdata.push_back(tmpstream.str());
                ++hh;
            }
            
            std::stringstream tmpstream;
            tmpstream << "currently_allocated: " ;
            for(unsigned int ii = 0;ii < INT_ATT_NAMES::enumArity();++ii)
            {
                std::string res = falsestring; 
                if (_currallocatedboatt[ii])
                    res = truestring;
                tmpstream << "att[" << ii << "]=" << res << " ";
            }
            _loginfo._currentlyallocated = tmpstream.str();
        }

        class EdgeVertInd
        {
        public:

            GLuint  _v[2];  // the two Vertex indices are ordered!
         
            EdgeVertInd() {}
            EdgeVertInd(const MESH_TYPE& m,typename MESH_TYPE::FacePointer  pf, const int nz) { this->set(m,pf,nz); }
            EdgeVertInd(const MESH_TYPE& m,typename MESH_TYPE::EdgePointer  pe, const int nz) { this->set(m,pe,nz); }

            void set(const MESH_TYPE& m,typename MESH_TYPE::FacePointer  pf, const int nz )
            {
                assert(pf!=0);
                assert(nz>=0);
                assert(nz<pf->VN());

                _v[0] = GLuint(vcg::tri::Index(m,pf->V(nz)));;
                _v[1] = GLuint(vcg::tri::Index(m,pf->V(pf->Next(nz))));
                assert(_v[0] != _v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

                if( _v[0] > _v[1] ) 
                    std::swap(_v[0],_v[1]);   
            }

            void set(const MESH_TYPE& m,typename MESH_TYPE::EdgePointer pe,const int nz)
            {
                assert(pe!=0);
                assert(nz>=0);
                assert(nz<2);

                _v[0] = size_t(vcg::tri::Index(m,pe->V(nz)));;
                _v[1] = size_t(vcg::tri::Index(m,pe->V((nz + 1)%2)));
                assert(_v[0] != _v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

                if( _v[0] > _v[1] )
                    std::swap(_v[0],_v[1]);   
            }

            inline bool operator<(const EdgeVertInd& pe) const
            {
                if(_v[0]<pe._v[0] ) 
                    return true;
                else if(_v[0]>pe._v[0] ) 
                    return false;
                else 
                    return _v[1] < pe._v[1];
            }

            inline bool operator==( const EdgeVertInd & pe ) const
            {
                return _v[0]==pe._v[0] && _v[1]==pe._v[1];
            }
        };

        static void fillEdgeVector(MESH_TYPE &m, std::vector<EdgeVertInd> &edgeVec, bool includeFauxEdge=true)
        {
            if (m.FN() > 0)
            {
                edgeVec.reserve(m.FN()*3);
                for(typename MESH_TYPE::FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
                    if( ! (*fi).IsD() )
                        for(int j=0;j<(*fi).VN();++j)
                            if(includeFauxEdge || !(*fi).IsF(j))
                                edgeVec.push_back(EdgeVertInd(m,&*fi,j));
            }
            else
                if ((m.VN() > 0) && (m.EN() > 0) )
                {
                    edgeVec.reserve(m.EN()*2);
                    for(typename MESH_TYPE::EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
                        if( ! (*ei).IsD() )
                            for(int j=0;j<2;++j)
                                    edgeVec.push_back(EdgeVertInd(m,&*ei,j));
                }

        }

        static void fillUniqueEdgeVector(MESH_TYPE &m, std::vector<EdgeVertInd> &edgeVec)
        {
            fillEdgeVector(m,edgeVec,false);
            std::sort(edgeVec.begin(), edgeVec.end());		

            typename std::vector<EdgeVertInd>::iterator newEnd = std::unique(edgeVec.begin(), edgeVec.end());

            edgeVec.resize(newEnd-edgeVec.begin());
        }

        struct GLBufferObject
        {
            GLBufferObject(size_t components,GLenum gltype,GLenum clientstatetag,GLenum target)
                :_size(0),_components(components),_isvalid(false),_gltype(gltype),_clientstatetag(clientstatetag),_target(target),_bohandle(0)
            {
            }

            GLBufferObject(size_t components,GLenum gltype,GLenum target)
                :_size(0),_components(components),_isvalid(false),_gltype(gltype),_clientstatetag(),_target(target),_bohandle(0)
            {
            }

            size_t getSizeOfGLType() const
            {
                switch(_gltype)
                {
                case(GL_FLOAT):
                    return sizeof(GLfloat);
                case(GL_INT):
                    return sizeof(GLint);
                case(GL_UNSIGNED_INT):
                    return sizeof(GLuint);
                case(GL_UNSIGNED_BYTE):
                    return sizeof(GLubyte);
                }
                return 0;
            }

            size_t _size;
            const size_t _components;
            bool _isvalid;
            const GLenum _gltype;
            const GLenum _target;

            /*WARNING!!!!!!!!!!!!!!!!! In openGL INDEX BO doesn't require to be enabled/disabled so has NOT a valid tag associated.
            In this case the client state tag remains not initialized and it's not meaningful */
            const GLenum _clientstatetag;
            /**********************************************************************************/

            GLuint _bohandle;
        };

        //ideally this should be const. I'm not yet sure if VCGLib will allow me to declare it as constant
        MESH_TYPE& _mesh;

        MemoryInfo& _gpumeminfo;

        /*The buffer objects used for the rendering operation. They are shared among all the views*/
        std::vector<GLBufferObject*> _bo;

        typedef std::map< UNIQUE_VIEW_ID_TYPE,PVData > ViewsMap;

        ///*_perviewreqatts contains a map of the requested atts by each single view. it's maintained for the actual rendering step*/
        ViewsMap _perviewreqatts;

        /*_currboatt contains the union of all the requested attributes by each single view on the scene. At the end it represents the BOs allocated in the GPU memory*/
        /* WARNING!!!! The currently allocated BOs are the union of all the BOs requested to be visualized in the _perviewreqatts plus, possibly, the edgeindex bo (depending by the kind of per view render primitive modality that is requested) and the vertexindex bo (depending of the solid rendering modality, per-face/per-vertex)*/
        /*             The EdgeIndex bo is allocated only if one of the requested rendering modality is PR_WIREFRAME_EDGES or PR_WIREFRAME_EDGES. If we have PR_SOLID the glPolygonMode function is used for rendering the triangle wireframe view*/ 
        InternalRendAtts _currallocatedboatt;

        bool _borendering;
        size_t _perbatchprim;

        /*Additional structures used for per wedge texturing modality*/ 
        typedef std::vector< std::pair< GLuint,GLuint > > ChunkVector;
        typedef std::map< short, ChunkVector > ChunkMap;

        std::vector< std::pair<short,GLuint> > _texindnumtriangles;
        ChunkMap  _chunkmap;
        
        //Horrible waste of memory space...but computing the list of edges is too much expensive...we must minimize it!
        std::vector<EdgeVertInd> _edge;
        size_t _meshverticeswhenedgeindiceswerecomputed;
        size_t _meshtriangleswhenedgeindiceswerecomputed;

        //vcg::GLOptions _glopts;
        vcg::Matrix44<typename MESH_TYPE::ScalarType> _tr;

        bool _debugmode;
        DebugInfo _loginfo;

        std::vector<InternalRendAtts> _meaningfulattsperprimitive;
    };
}

#endif
