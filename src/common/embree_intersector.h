#pragma once


#include <Eigen/Geometry>
#include <Eigen/Core>

#include <embree4/rtcore.h>
#include <embree4/rtcore_ray.h>
#include <iostream>
#include <vector>

// Reimplementation of the embree::Hit struct from embree1.0
//
// TODO: template on floating point type
struct Hit
{
    int id; // primitive id
    int gid; // geometry id (not used)
    // barycentric coordinates so that
    //   pos = V.row(F(id,0))*(1-u-v)+V.row(F(id,1))*u+V.row(F(id,2))*v;
    float u,v;
    // parametric distance so that
    //   pos = origin + t * dir
    float t;
};


// keep track of embree device
struct EmbreeDevice
{
    RTCDevice embree_device;
    int       embree_device_cntr;

    static EmbreeDevice & instance()
    {
        static EmbreeDevice s;
        return s;
    } // instance

    EmbreeDevice(const EmbreeDevice &) = delete;
    EmbreeDevice & operator = (const EmbreeDevice &) = delete;

    static RTCDevice get_device(const char *config=nullptr)
    {
        return instance().get(config);
    }

    static void release_device(void)
    {
        instance().release();
    }

private:

    EmbreeDevice():embree_device(nullptr),embree_device_cntr(0) {}

    ~EmbreeDevice()
    {
        if(embree_device)
            rtcReleaseDevice(embree_device);
    }

    RTCDevice
    get(const char *config=nullptr)
    {
        if(!embree_device)
        {
            embree_device = rtcNewDevice (config);
            if(rtcGetDeviceError (embree_device) != RTC_ERROR_NONE)
                std::cerr << "Embree: An error occurred while initializing embree core!" << std::endl;
        }
        ++embree_device_cntr;
        return embree_device;
    }

    void release()
    {
        if(!--embree_device_cntr) {
            rtcReleaseDevice (embree_device);
            embree_device = nullptr;
        }
    }
};




class EmbreeIntersector
{
public:
    typedef Eigen::Matrix<float,Eigen::Dynamic,3> PointMatrixType;
    typedef Eigen::Matrix<int,Eigen::Dynamic,3> FaceMatrixType;
public:
    EmbreeIntersector();
private:
    // Copying and assignment are not allowed.
    EmbreeIntersector(const EmbreeIntersector & that);
    EmbreeIntersector & operator=(const EmbreeIntersector &);
public:
    virtual ~EmbreeIntersector();

    // Initialize with a given mesh.
    //
    // Inputs:
    //   V  #V by 3 list of vertex positions
    //   F  #F by 3 list of Oriented triangles
    //   isStatic  scene is optimized for static geometry
    // Side effects:
    //   The first time this is ever called the embree engine is initialized.
    void init(
            const PointMatrixType& V,
            const FaceMatrixType& F,
            bool isStatic = false);

    // Initialize with a given mesh.
    //
    // Inputs:
    //   V  vector of #V by 3 list of vertex positions for each geometry
    //   F  vector of #F by 3 list of Oriented triangles for each geometry
    //   masks  a 32 bit mask to identify active geometries.
    //   isStatic  scene is optimized for static geometry
    // Side effects:
    //   The first time this is ever called the embree engine is initialized.
    void init(
            const std::vector<const PointMatrixType*>& V,
            const std::vector<const FaceMatrixType*>& F,
            const std::vector<int>& masks,
            bool isStatic = false);

    // Deinitialize embree datasctructures for current mesh.  Also called on
    // destruction: no need to call if you just want to init() once and
    // destroy.
    void deinit();

    // Given a ray find the first hit
    //
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    //   tnear      start of ray segment
    //   tfar       end of ray segment
    //   masks      a 32 bit mask to identify active geometries.
    // Output:
    //   hit        information about hit
    // Returns true if and only if there was a hit
    bool intersectRay(
            const Eigen::RowVector3f& origin,
            const Eigen::RowVector3f& direction,
            Hit& hit,
            float tnear = 0,
            float tfar = std::numeric_limits<float>::infinity(),
            int mask = 0xFFFFFFFF) const;

    // Given a ray find the first hit
    // This is a conservative hit test where multiple rays within a small radius
    // will be tested and only the closesest hit is returned.
    //
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    //   tnear      start of ray segment
    //   tfar       end of ray segment
    //   masks      a 32 bit mask to identify active geometries.
    //   geoId      id of geometry mask (default std::numeric_limits<float>::infinity() if no: no masking)
    //   closestHit true for gets closest hit, false for furthest hit
    // Output:
    //   hit        information about hit
    // Returns true if and only if there was a hit
    bool intersectBeam(
            const Eigen::RowVector3f& origin,
            const Eigen::RowVector3f& direction,
            Hit& hit,
            float tnear = 0,
            float tfar = std::numeric_limits<float>::infinity(),
            int mask = 0xFFFFFFFF,
            int geoId = -1,
            bool closestHit = true,
            unsigned int samples = 4) const;

    // Given a ray find all hits in order
    //
    // Inputs:
    //   origin     3d origin point of ray
    //   direction  3d (not necessarily normalized) direction vector of ray
    //   tnear      start of ray segment
    //   tfar       end of ray segment
    //   masks      a 32 bit mask to identify active geometries.
    // Output:
    //   hit        information about hit
    //   num_rays   number of rays shot (at least one)
    // Returns true if and only if there was a hit
    bool intersectRay(
            const Eigen::RowVector3f& origin,
            const Eigen::RowVector3f& direction,
            std::vector<Hit > &hits,
            int& num_rays,
            float tnear = 0,
            float tfar = std::numeric_limits<float>::infinity(),
            int mask = 0xFFFFFFFF) const;

    // Given a ray find the first hit
    //
    // Inputs:
    //   a    3d first end point of segment
    //   ab   3d vector from a to other endpoint b
    // Output:
    //   hit  information about hit
    // Returns true if and only if there was a hit
    bool intersectSegment(
            const Eigen::RowVector3f& a,
            const Eigen::RowVector3f& ab,
            Hit &hit,
            int mask = 0xFFFFFFFF) const;

private:

    struct Vertex   {float x,y,z,a;};
    struct Triangle {int v0, v1, v2;};

    RTCScene scene;
    unsigned geomID;
    Vertex* vertices;
    Triangle* triangles;
    bool initialized;

    RTCDevice device;

    void createRay(
            RTCRayHit& ray,
            const Eigen::RowVector3f& origin,
            const Eigen::RowVector3f& direction,
            float tnear,
            float tfar,
            int mask) const;
};

