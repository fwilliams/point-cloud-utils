#ifndef GL_FIELD
#define GL_FIELD

#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <vcg/complex/allocate.h>

namespace vcg{
template <class MeshType>
class GLField
{
	typedef typename MeshType::FaceType FaceType;
	typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
	typedef typename MeshType::ScalarType ScalarType;
	
public:

	static void GLDrawField(CoordType dir[4],
							CoordType center,
                            ScalarType &size,
                            bool onlyPD1=false,
                            bool oneside=false)
	{
        ScalarType size1=size;
        ScalarType size2=size;
        if (oneside)size2=0;

        glLineWidth(2);
        //vcg::glColor(vcg::Color4b(0,0,255,255));
        vcg::glColor(vcg::Color4b(0,0,0,255));
        glBegin(GL_LINES);
            glVertex(center+dir[0]*size1);
            glVertex(center-dir[0]*size2);
        glEnd();

        if (onlyPD1)return;
        glLineWidth(2);
        //vcg::glColor(vcg::Color4b(0,255,0,255));
        vcg::glColor(vcg::Color4b(0,0,0,255));
        glBegin(GL_LINES);
            glVertex(center+dir[1]*size1);
            glVertex(center-dir[1]*size2);
        glEnd();

	}

    ///draw the cross field of a given face in a given position
    static void GLDrawSingleFaceField(const FaceType &f,
                                CoordType pos,
                                ScalarType &size,
                                bool onlyPD1=false,
                                bool oneside=false)
    {
        CoordType center=pos;
        CoordType normal=f.cN();
        CoordType dir[4];
        vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
        GLDrawField(dir,center,size,onlyPD1,oneside);
    }

	///draw the cross field of a given face
    static void GLDrawSingleFaceField(const FaceType &f,
                                ScalarType &size,
                                bool onlyPD1=false,
                                bool oneside=false)
	{
        CoordType center=(f.cP(0)+f.cP(1)+f.cP(2))/3;
		CoordType normal=f.cN();
		CoordType dir[4];
		vcg::tri::CrossField<MeshType>::CrossVector(f,dir);
        GLDrawField(dir,center,size,onlyPD1,oneside);
	}
	
//    static void GLDrawFaceSeams(const FaceType &f,
//                                vcg::Point3<bool> seams,
//                                vcg::Color4b seamCol[3])
//    {
//        glLineWidth(2);

//        glBegin(GL_LINES);
//        for (int i=0;i<3;i++)
//        {
//            if (!seams[i])continue;
//            vcg::glColor(seamCol[i]);
//            glVertex(f.V0(i)->P());
//            glVertex(f.V1(i)->P());
//        }
//        glEnd();
//    }

    static void GLDrawVertField(const VertexType &v,
                                ScalarType &size)
    {
        CoordType center=v.cP();
        CoordType normal=v.cN();
        CoordType dir[4];
        vcg::tri::CrossField<MeshType>::CrossVector(v,dir);
        GLDrawField(dir,center,size);
    }


    static void GLDrawFaceField(const MeshType &mesh,
                                bool onlyPD1=false,
                                bool oneside=false,
                                ScalarType scale=0.002)
	{

		glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDepthRange(0.0,0.999);
		glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);
        ScalarType size=mesh.bbox.Diag()*scale;
        for (unsigned int i=0;i<mesh.face.size();i++)
		{
            if (mesh.face[i].IsD())continue;
            GLDrawSingleFaceField(mesh.face[i],size,onlyPD1,oneside);
		}
		glPopAttrib();
	}

    static void GLDrawVertField(const MeshType &mesh,ScalarType sizeF=0.01)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDepthRange(0.0,0.9999);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);
        ScalarType size=mesh.bbox.Diag()*sizeF;
        for (int i=0;i<mesh.vert.size();i++)
		{
            if (mesh.vert[i].IsD())continue;
            GLDrawVertField(mesh.vert[i],size);
		}
		glPopAttrib();
	}

    static void GLDrawSingularity(MeshType &mesh)
    {
        // query if an attribute is present or not
       bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
       bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularIndex"));

       if (!hasSingular)return;
       if(!hasSingularIndex)return;

       typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
       Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));
       typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;
       Handle_SingularIndex =vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularIndex"));

       glPushAttrib(GL_ALL_ATTRIB_BITS);

       glDepthRange(0.0,0.9999);
       glEnable(GL_COLOR_MATERIAL);
       glDisable(GL_LIGHTING);
       glDisable(GL_BLEND);
       glPointSize(20);
       glBegin(GL_POINTS);
       for (size_t i=0;i<mesh.vert.size();i++)
       {
           if (mesh.vert[i].IsD())continue;
           if (!Handle_Singular[i])continue;


           int SingIndex=Handle_SingularIndex[i];

           vcg::Color4b colSing;

           switch (SingIndex)
           {
             case 1:colSing=vcg::Color4b(0,0,255,255);      break;
             case 2:colSing=vcg::Color4b(0,255,0,255);    break;
             case 3:colSing=vcg::Color4b(255,0,0,255);      break;
             case 4:colSing=vcg::Color4b(255,255,0,255);      break;
             default:colSing=vcg::Color4b(255,0,255,255);
           }


           vcg::glColor(colSing);
           vcg::glVertex(mesh.vert[i].P());
       }
       glEnd();
       glPopAttrib();
    }
};

}

#endif
