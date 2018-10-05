#include <iostream>
#include <QTime>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "nanoflann.hpp"

#include <vcg/complex/complex.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>


#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/perfect_spatial_hashing.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/space/index/octree.h>

int num_test = 1000;
int kNearest = 256;
float queryDist = 0.0037;
float ratio = 1000.0f;


class CVertex;
class CFace;
class CEdge;

class CUsedTypes	: public vcg::UsedTypes < vcg::Use< CVertex >::AsVertexType, vcg::Use< CFace >::AsFaceType>{};
class CVertex		: public vcg::Vertex < CUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Radiusf, vcg::vertex::BitFlags, vcg::vertex::Qualityf, vcg::vertex::Color4b>{};
class CFace			: public vcg::Face < CUsedTypes, vcg::face::VertexRef>{};

class CMesh			: public vcg::tri::TriMesh < std::vector< CVertex >, std::vector< CFace > > {};


template <typename T>
struct PointCloud
{
    struct Point
    {
        T  x,y,z;
    };

    std::vector<Point>  pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t size) const
    {
        const T d0=p1[0]-pts[idx_p2].x;
        const T d1=p1[1]-pts[idx_p2].y;
        const T d2=p1[2]-pts[idx_p2].z;
        return d0*d0+d1*d1+d2*d2;
    }

    inline T kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0) return pts[idx].x;
        else if (dim==1) return pts[idx].y;
        else return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }

};


void testKDTree(CMesh& mesh, std::vector<unsigned int>& test_indeces, std::vector<vcg::Point3f>& randomSamples)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "KDTree" << std::endl;
    QTime time;
    time.start();

    // Construction of the kdTree
    vcg::ConstDataWrapper<CMesh::VertexType::CoordType> wrapperVcg(&mesh.vert[0].P(), mesh.vert.size(), size_t(mesh.vert[1].P().V()) - size_t(mesh.vert[0].P().V()));
    vcg::KdTree<CMesh::ScalarType> kdTreeVcg(wrapperVcg);
    std::cout << "Build: " << time.elapsed() <<  " ms" << std::endl;
    int nn=1;
    // Computation of the point radius
    float mAveragePointSpacing = 0;
    time.restart();
    #pragma omp parallel for reduction(+: mAveragePointSpacing) schedule(dynamic, 10)
    for (int i = 0; i < mesh.vert.size(); i++)
    {
#ifdef #ifdef _OPENMP
      nn =omp_get_num_threads();
#endif
        vcg::KdTree<CMesh::ScalarType>::PriorityQueue queue;
        kdTreeVcg.doQueryK(mesh.vert[i].cP(), 16, queue);
        float newRadius = 2.0f * sqrt(queue.getWeight(0)/ queue.getNofElements());
        mesh.vert[i].R() -= newRadius;
        mAveragePointSpacing += newRadius;
    }
    std::cout << "Num trhread " << nn << std::endl;
    mAveragePointSpacing /= mesh.vert.size();
    std::cout << "Average point radius (OpenMP with" << nn << " threads) " << mAveragePointSpacing << std::endl;
    std::cout << "Time (OpenMP): " << time.elapsed() << " ms" << std::endl;

    queryDist = mAveragePointSpacing * 150;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float avgTime = 0.0f;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        std::vector<unsigned int> indeces;
        std::vector<float> dists;
        kdTreeVcg.doQueryDist(mesh.vert[test_indeces[ii]].cP(), queryDist, indeces, dists);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl;

    // Test with the k-nearest search
    std::cout << "k-Nearest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        vcg::KdTree<CMesh::ScalarType>::PriorityQueue queue;
        kdTreeVcg.doQueryK(mesh.vert[test_indeces[ii]].cP(), kNearest, queue);
        avgTime += time.elapsed();
    }
    std::cout << "Time (k = " << kNearest << "): " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl;

    // Test with the closest search
    std::cout << "Closest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        unsigned int index;
        float minDist;
        kdTreeVcg.doQueryClosest(randomSamples[ii], index, minDist);
        avgTime += time.elapsed();
    }
    std::cout << "Time : " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl << std::endl;
}


void testNanoFLANN(CMesh& mesh, std::vector<unsigned int>& test_indeces, std::vector<vcg::Point3f> randomSamples)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "nanoFLANN" << std::endl;

    PointCloud<float> cloud;
    cloud.pts.resize(mesh.vert.size());
    for (size_t i=0; i < mesh.vert.size(); i++)
    {
        cloud.pts[i].x = mesh.vert[i].P().X();
        cloud.pts[i].y = mesh.vert[i].P().Y();
        cloud.pts[i].z = mesh.vert[i].P().Z();
    }

    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, PointCloud<float> > ,
        PointCloud<float>,
        3 /* dim */
        > my_kd_tree_t;

    // Construction of the nanoFLANN KDtree
    QTime time;
    time.start();
    my_kd_tree_t   index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(16) );
    index.buildIndex();
    std::cout << "Build nanoFlann: " << time.elapsed() <<  " ms" << std::endl;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float avgTime = 0.0f;
    std::vector<std::pair<size_t,float> >   ret_matches;
    nanoflann::SearchParams params;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        const size_t nMatches = index.radiusSearch(mesh.vert[test_indeces[ii]].P().V(), queryDist, ret_matches, params);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl;

    // Test with the k-nearest search
    std::cout << "k-Nearest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    std::vector<size_t>   ret_index(kNearest);
    std::vector<float> out_dist_sqr(kNearest);
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        index.knnSearch(mesh.vert[test_indeces[ii]].P().V(), kNearest, &ret_index[0], &out_dist_sqr[0]);
        avgTime += time.elapsed();
    }
    std::cout << "Time (k = " << kNearest << "): " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl;

    // Test with the closest search
    std::cout << "Closest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    std::vector<size_t>   ret_index_clos(1);
    std::vector<float> out_dist_sqr_clos(1);
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        index.knnSearch(randomSamples[ii].V(), 1, &ret_index_clos[0], &out_dist_sqr_clos[0]);
        avgTime += time.elapsed();
    }
    std::cout << "Time : " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl << std::endl;
}


void testUniformGrid(CMesh& mesh, std::vector<unsigned int>& test_indeces, std::vector<vcg::Point3f>& randomSamples)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "Uniform Grid" << std::endl;
    QTime time;
    time.start();

    // Construction of the uniform grid
    typedef vcg::GridStaticPtr<CMesh::VertexType, CMesh::VertexType::ScalarType> MeshGrid;
    MeshGrid uniformGrid;
    uniformGrid.Set(mesh.vert.begin(), mesh.vert.end());
    std::cout << "Build: " << time.elapsed() <<  " ms" << std::endl;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float  avgTime = 0.0f;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetInSphereVertex(mesh, uniformGrid, mesh.vert[test_indeces[ii]].cP(), queryDist, vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl;

    // Test with the k-nearest search
    std::cout << "k-Nearest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetKClosestVertex(mesh, uniformGrid, kNearest, mesh.vert[test_indeces[ii]].cP(), mesh.bbox.Diag(), vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (k = " << kNearest << "): " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl;

    // Test with the Closest search
    std::cout << "Closest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        float minDist;
        vcg::tri::GetClosestVertex(mesh, uniformGrid, randomSamples[ii], mesh.bbox.Diag(), minDist);
        avgTime += time.elapsed();
    }
    std::cout << "Time : " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl << std::endl;
}



void testSpatialHashing(CMesh& mesh, std::vector<unsigned int>& test_indeces, std::vector<vcg::Point3f>& randomSamples)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "Spatial Hashing" << std::endl;
    QTime time;
    time.start();

    // Construction of the uniform grid
    typedef vcg::SpatialHashTable<CMesh::VertexType, CMesh::VertexType::ScalarType> MeshGrid;
    MeshGrid uniformGrid;
    uniformGrid.Set(mesh.vert.begin(), mesh.vert.end());
    std::cout << "Build: " << time.elapsed() <<  " ms" << std::endl;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float  avgTime = 0.0f;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetInSphereVertex(mesh, uniformGrid, mesh.vert[test_indeces[ii]].cP(), queryDist, vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl;

    // Test with the k-nearest search
    std::cout << "k-Nearest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetKClosestVertex(mesh, uniformGrid, kNearest, mesh.vert[test_indeces[ii]].cP(), mesh.bbox.Diag(), vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (k = " << kNearest << "): " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl;

    // Test with the Closest search
    std::cout << "Closest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        float minDist;
        vcg::tri::GetClosestVertex(mesh, uniformGrid, randomSamples[ii], mesh.bbox.Diag(), minDist);
        avgTime += time.elapsed();
    }
    std::cout << "Time : " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl << std::endl;
}



void testPerfectSpatialHashing(CMesh& mesh, std::vector<unsigned int>& test_indeces)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "Perfect Spatial Hashing" << std::endl;
    QTime time;
    time.start();

    // Construction of the uniform grid
    typedef vcg::SpatialHashTable<CMesh::VertexType, CMesh::VertexType::ScalarType> MeshGrid;
    MeshGrid uniformGrid;
    uniformGrid.Set(mesh.vert.begin(), mesh.vert.end());
    std::cout << "Build: " << time.elapsed() <<  " ms" << std::endl;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float  avgTime = 0.0f;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetInSphereVertex(mesh, uniformGrid, mesh.vert[test_indeces[ii]].cP(), queryDist, vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl << std::endl;
}


void testOctree(CMesh& mesh, std::vector<unsigned int>& test_indeces, std::vector<vcg::Point3f>& randomSamples)
{
    std::cout << "==================================================="<< std::endl;
    std::cout << "Octree" << std::endl;
    QTime time;
    time.start();

    // Construction of the uniform grid
    typedef vcg::Octree<CMesh::VertexType, CMesh::VertexType::ScalarType> MeshGrid;
    MeshGrid uniformGrid;
    uniformGrid.Set(mesh.vert.begin(), mesh.vert.end());
    std::cout << "Build: " << time.elapsed() <<  " ms" << std::endl;

    // Test with the radius search
    std::cout << "Radius search (" << num_test << " tests)"<< std::endl;
    float  avgTime = 0.0f;
    for (int ii = 0; ii < num_test; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetInSphereVertex(mesh, uniformGrid, mesh.vert[test_indeces[ii]].cP(), queryDist, vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (radius = " << queryDist << "): " << avgTime << " ms (mean " << avgTime / num_test << "ms)"  << std::endl;

    // Test with the k-nearest search
    std::cout << "k-Nearest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        std::vector<CMesh::VertexPointer> vertexPtr;
        std::vector<CMesh::VertexType::CoordType> points;
        std::vector<float> dists;
        vcg::tri::GetKClosestVertex(mesh, uniformGrid, kNearest, mesh.vert[test_indeces[ii]].cP(), mesh.bbox.Diag(), vertexPtr, dists, points);
        avgTime += time.elapsed();
    }
    std::cout << "Time (k = " << kNearest << "): " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl;

    // Test with the Closest search
    std::cout << "Closest search (" << num_test*10 << " tests)"<< std::endl;
    avgTime = 0.0f;
    for (int ii = 0; ii < num_test * 10; ii++)
    {
        time.restart();
        float minDist;
        vcg::tri::GetClosestVertex(mesh, uniformGrid, randomSamples[ii], mesh.bbox.Diag(), minDist);
        avgTime += time.elapsed();
    }
    std::cout << "Time : " << avgTime << " ms (mean " << avgTime / (num_test * 10) << "ms)"  << std::endl << std::endl;
}



int main( int argc, char * argv[] )
{
    if (argc < 2) {
        std::cout << "Invalid arguments" << std::endl;
        exit(-1);
    }
    CMesh mesh;
    if (vcg::tri::io::Importer<CMesh>::Open(mesh, argv[1])  != 0)
        std::cout << "Invalid filename" << std::endl;

    std::cout << "Mesh BBox diagonal: " << mesh.bbox.Diag() << std::endl;
    std::cout << "Max point random offset: " << mesh.bbox.Diag() / 1000.0f << std::endl << std::endl;

    vcg::math::MarsenneTwisterRNG randGen;
    randGen.initialize(0);
    std::vector<vcg::Point3f> randomSamples;
    for (int i = 0; i < num_test * 10; i++)
        randomSamples.push_back(vcg::math::GeneratePointOnUnitSphereUniform<float>(randGen) * randGen.generate01() * mesh.bbox.Diag() / ratio);

    std::vector<unsigned int> test_indeces;
    for (int i = 0; i < num_test * 10; i++)
    {
        int index = randGen.generate01() * (mesh.vert.size() - 1);
        test_indeces.push_back(index);
        randomSamples[i] += mesh.vert[i].P();
    }

    testKDTree(mesh, test_indeces, randomSamples);
    testNanoFLANN(mesh, test_indeces, randomSamples);
    testUniformGrid(mesh, test_indeces, randomSamples);
    testSpatialHashing(mesh, test_indeces, randomSamples);
    testPerfectSpatialHashing(mesh, test_indeces);
    testOctree(mesh, test_indeces, randomSamples);
}
