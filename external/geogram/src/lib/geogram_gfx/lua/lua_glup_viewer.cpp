/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifdef __ANDROID__

#include <geogram/lua/lua_wrap.h>

void init_lua_glup_viewer(lua_State* L) {
    (void)L;
}

#else

#include <geogram_gfx/lua/lua_glup_viewer.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/lua/lua_wrap.h>
#include <geogram/basic/stopwatch.h>

namespace {
    using namespace GEO;

    namespace LUAGLUPVIEWERImpl {
	static double t0 = 0.0;
	
	static int ElapsedTime(lua_State* L) {
	    if(lua_gettop(L) != 0) {
		return luaL_error(
		    L, "'GLUP.ElapsedTime()' invalid number of arguments"
		);
	    }
	    double result = 0.0;
	    result = GEO::SystemStopwatch::now() - t0;
	    lua_pushnumber(L,double(result));    
	    return 1;
	}

	static int ResetViewer(lua_State* L) {
	    if(lua_gettop(L) != 0) {
		return luaL_error(
		    L, "'GLUP.ResetViewer()' invalid number of arguments"
		);
	    }
	    glup_viewer_home();
	    glup_viewer_disable(GLUP_VIEWER_CLIP);
	    glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
	    glup_viewer_disable(GLUP_VIEWER_FIXED_CLIP);    
	    glup_viewer_disable(GLUP_VIEWER_ROTATE_LIGHT);
	    glup_viewer_enable(GLUP_VIEWER_BACKGROUND);    
	    GEO::Application::instance()->set_lighting(true);
	    GEO::Application::instance()->set_background_color_1(
		1.0f, 1.0f, 1.0f
	    );
	    GEO::Application::instance()->set_background_color_2(
		0.0f, 0.0f, 0.7f
	    );	    
	    t0 = GEO::SystemStopwatch::now();
	    return 0;
	}

	static int ArcadeStyle(lua_State* L) {
	    if(lua_gettop(L) != 0) {
		return luaL_error(
		    L, "'GLUP.ArcadeStyle()' invalid number of arguments"
		);
	    }
	    glup_viewer_home();
	    glup_viewer_disable(GLUP_VIEWER_CLIP);
	    glup_viewer_disable(GLUP_VIEWER_EDIT_CLIP);
	    glup_viewer_disable(GLUP_VIEWER_FIXED_CLIP);    
	    glup_viewer_disable(GLUP_VIEWER_ROTATE_LIGHT);
	    glup_viewer_disable(GLUP_VIEWER_BACKGROUND);    
	    GEO::Application::instance()->set_lighting(false);
	    GEO::Application::instance()->set_background_color_1(0,0,0);
	    GEO::Application::instance()->set_background_color_2(0,0,0);
	    return 0;    
	}

	static int SetRegionOfInterest(lua_State* L) {
	    if(lua_gettop(L) != 6) {
		return luaL_error(
		    L,
		    "'GLUP.SetRegionOfInterest()' invalid number of arguments"
		);
	    }
	    if(
		!lua_isnumber(L,1) ||
		!lua_isnumber(L,2) ||
		!lua_isnumber(L,3) ||
		!lua_isnumber(L,4) ||	
		!lua_isnumber(L,5) ||
		!lua_isnumber(L,6) 
	    ) {
		return luaL_error(
		    L, "'GLUP.SetRegionOfInterest()' arguments should be numbers"
		);
	    }
	    glup_viewer_set_region_of_interest(
		float(lua_tonumber(L,1)),
		float(lua_tonumber(L,2)),
		float(lua_tonumber(L,3)),
		float(lua_tonumber(L,4)),
		float(lua_tonumber(L,5)),
		float(lua_tonumber(L,6))	
	    );
	    return 0;
	}

	static int GetRegionOfInterest(lua_State* L) {
	    if(lua_gettop(L) != 0) {
		return luaL_error(
		    L, "'GLUP.GetRegionOfInterest()' invalid number of arguments"
		);
	    }
	    float xm,ym,zm,xM,yM,zM;
	    glup_viewer_get_region_of_interest(
		&xm, &ym, &zm, &xM, &yM, &zM
	    );
	    lua_pushnumber(L,double(xm));
	    lua_pushnumber(L,double(ym));
	    lua_pushnumber(L,double(zm));
	    lua_pushnumber(L,double(xM));
	    lua_pushnumber(L,double(yM));
	    lua_pushnumber(L,double(zM));            
	    return 6;
	}
    }
}
    
void init_lua_glup_viewer(lua_State* L) {
    lua_getglobal(L,"GLUP");
    geo_assert(!lua_isnil(L,-1)); // Make sure GLUP was registered before.
    GEO::lua_bindwrapper(L,LUAGLUPVIEWERImpl::ElapsedTime);
    GEO::lua_bindwrapper(L,LUAGLUPVIEWERImpl::ResetViewer);
    GEO::lua_bindwrapper(L,LUAGLUPVIEWERImpl::ArcadeStyle);
    GEO::lua_bindwrapper(L,LUAGLUPVIEWERImpl::SetRegionOfInterest);
    GEO::lua_bindwrapper(L,LUAGLUPVIEWERImpl::GetRegionOfInterest);
    lua_pop(L,1);
}

#endif
