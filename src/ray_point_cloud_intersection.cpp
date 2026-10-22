#include <npe.h>
#include <common/embree_intersector.h>
#include <tuple>
#include <numeric>

#include "common/common.h"
#include "common/strutil.h"


namespace {
    enum class GeometryType {
        SPHERE,
        CIRCLE,
    };

    GeometryType geometry_type_from_string(const std::string& geometry_type_str) {
        if (strutil::compare_ignore_case(geometry_type_str, "sphere")) {
            return GeometryType::SPHERE;
        } else if (strutil::compare_ignore_case(geometry_type_str, "circle")) {
            return GeometryType::CIRCLE;
        } else {
            throw pybind11::value_error("Invalid geometry_type. Got " + geometry_type_str +
                                        ". But expected 'sphere' or 'circle'");
        }
    }

    template <typename T1, typename T2>
    bool validate_rays(T1 ray_o, T2 ray_d) {
        bool use_single_ray_origin = false;
        if (ray_o.size() == 3) {
            use_single_ray_origin = true;
        } else {
            if (ray_o.rows() != ray_d.rows()) {
                throw pybind11::value_error("ray_o and ray_d must have the same number of rows (one ray origin per ray direction). "
                                            "(Note: ray_o can have one row to use the same origin for all directions)");
            }
        }
        if (ray_o.cols() != 3 && !use_single_ray_origin) {
            throw pybind11::value_error("Invalid shape for ray_o must have shape (N, 3) but got (" +
                                        std::to_string(ray_o.rows()) + ", " + std::to_string(ray_o.cols()) + ").");
        }
        if (ray_d.cols() != 3) {
            throw pybind11::value_error("Invalid shape for ray_d must have shape (N, 3) but got (" +
                                        std::to_string(ray_d.rows()) + ", " + std::to_string(ray_d.cols()) + ").");
        }
        return use_single_ray_origin;
    }

    template <typename T1, typename T2, typename T3>
    GeometryType validate_point_geometry(T1 v, T2 n, T3 geometry_radius, int geometry_subdivisions_1, int geometry_subdivisions_2, std::string geometry_type) {
        validate_point_cloud_normals(v, n);

        if (geometry_radius.rows() != v.rows() || geometry_radius.cols() != 1) {
            throw pybind11::value_error("Invalid shape for geometry_radius, must have one row per vertex.");
        }

        if (geometry_subdivisions_1 < 4) {
            throw pybind11::value_error("Invalid geometry_subdivisions_1 is less than or equal to 4.");
        }

        if (geometry_subdivisions_2 < 4) {
            throw pybind11::value_error("Invalid geometry_subdivisions_1 is less than or equal to 4.");
        }

        return geometry_type_from_string(geometry_type);
    }

    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    using Matrix3 = Eigen::Matrix<T , 3, 3, Eigen::RowMajor>;

    template <typename FloatScalar>
    Matrix3<FloatScalar> local_basis(Vector3<FloatScalar> normal) {
        using V3 = Vector3<FloatScalar>;
        V3 ni = normal / normal.norm();
        V3 right;
        if (fabs(fabs(V3(0.0, 1.0, 0.0).dot(ni)) - 1.0) < 1e-5) {
            right = ni.cross(V3(1.0, 0.0, 0.0));
        } else {
            right = ni.cross(V3(0.0, 1.0, 0.0));
        }
        right = right / right.norm();
        V3 up = ni.cross(right);
        up = up / up.norm();
        Matrix3<FloatScalar> R;
        R << right[0], up[0], ni[0],
             right[1], up[1], ni[1],
             right[2], up[2], ni[2];
        return R;
    }


    template <typename OutScalar>
    void make_sphere_geometry(int stackCount, int sectorCount, double radius, int indexOffset, int vertexOffset,
                              double baseX, double baseY, double baseZ,
                              Eigen::Matrix<OutScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& vertices,
                              Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& indices) {
        double x, y, z, xy;                              // vertex position
        double nx, ny, nz, lengthInv = 1.0 / radius;    // vertex normal
        double s, t;                                     // vertex texCoord

        double sectorStep = 2 * M_PI / sectorCount;
        double stackStep = M_PI / stackCount;
        double sectorAngle, stackAngle;

        int currentVertex = vertexOffset;
        for(int i = 0; i <= stackCount; ++i) {
            stackAngle = M_PI / 2 - i * stackStep;      // starting from pi/2 to -pi/2
            xy = radius * cos(stackAngle);             // r * cos(u)
            z = radius * sin(stackAngle);              // r * sin(u)

            // add (sectorCount+1) vertices per stack
            // the first and last vertices have same position and normal, but different tex coords
            for(int j = 0; j <= sectorCount; ++j)
            {
                sectorAngle = j * sectorStep;           // starting from 0 to 2pi

                // vertex position (x, y, z)
                x = xy * cos(sectorAngle);             // r * cos(u) * cos(v)
                y = xy * sin(sectorAngle);             // r * cos(u) * sin(v)
                vertices(currentVertex, 0) = (OutScalar) baseX + x;
                vertices(currentVertex, 1) = (OutScalar) baseY + y;
                vertices(currentVertex, 2) = (OutScalar) baseZ + z;
                currentVertex += 1;
            }
        }

        int currentIndex = indexOffset;
        int k1, k2;
        for(int i = 0; i < stackCount; ++i)
        {
            k1 = i * (sectorCount + 1);     // beginning of current stack
            k2 = k1 + sectorCount + 1;      // beginning of next stack

            for(int j = 0; j < sectorCount; ++j, ++k1, ++k2)
            {
                // 2 triangles per sector excluding first and last stacks
                // k1 => k2 => k1+1
                if(i != 0)
                {
                    indices(currentIndex, 0) = vertexOffset + k1;
                    indices(currentIndex, 1) = vertexOffset + k2;
                    indices(currentIndex, 2) = vertexOffset + k1 + 1;
                    currentIndex += 1;
                }

                // k1+1 => k2 => k2+1
                if(i != (stackCount-1))
                {
                    indices(currentIndex, 0) = vertexOffset + k1 + 1;
                    indices(currentIndex, 1) = vertexOffset + k2;
                    indices(currentIndex, 2) = vertexOffset + k2 + 1;
                    currentIndex += 1;
                }
            }
        }
    }

    template <typename TypeV, typename TypeN, typename TypeR, typename OutScalar>
    int generate_splat_geometry(GeometryType geometry_type, int geometry_subdivisions_1, int geometry_subdivisions_2,
                                const TypeV& v, const TypeN& n, const TypeR& geometry_radius,
                                Eigen::Matrix<OutScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_v,
                                Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_f) {

        int num_vertices_per_geom = 0;
        int num_faces_per_geom = 0;
        switch (geometry_type) {
            case GeometryType::SPHERE: {
                num_vertices_per_geom = (geometry_subdivisions_1 + 1) * (geometry_subdivisions_2 + 1);
                num_faces_per_geom = (geometry_subdivisions_1 - 1) * geometry_subdivisions_2 * 2;
                break;
            }
            case GeometryType::CIRCLE: {
                num_vertices_per_geom = geometry_subdivisions_1 + 1;
                num_faces_per_geom = geometry_subdivisions_1;
                break;
            }
            default: {
                throw pybind11::value_error("Invalid geometry_type.");
                break;
            }
        }

        out_v.resize(v.rows() * num_vertices_per_geom, 3);
        out_f.resize(v.rows() * num_faces_per_geom, 3);

        using FloatScalar = OutScalar;
        using V3 = Vector3<FloatScalar>;
        using M3 = Matrix3<FloatScalar>;

        for (int i = 0; i < v.rows(); i += 1) {
            const int base_offset_v = i * num_vertices_per_geom;
            const int base_offset_f = i * num_faces_per_geom;
            switch (geometry_type) {
                case GeometryType::CIRCLE: {
                    const V3 ni((FloatScalar) n(i, 0), (FloatScalar) n(i, 1), (FloatScalar) n(i, 2));
                    const V3 vi((FloatScalar) v(i, 0), (FloatScalar) v(i, 1), (FloatScalar) v(i, 2));
                    const M3 R = local_basis(ni);

                    for (int j = 0; j < num_vertices_per_geom - 1; j += 1) {
                        const float theta = ((float) j / ((float) num_vertices_per_geom - 1.0f)) * (float) M_PI * 2.0f;
                        const float x0 = geometry_radius(i, 0) * cos(theta),
                                y0 = geometry_radius(i, 0) * sin(theta),
                                z0 = 0.0f;
                        const V3 vj = R * V3(x0, y0, 0.0) + vi;
                        for (int k = 0; k < 3; k += 1) {
                            out_v(base_offset_v + j, k) = vj[k];
                        }
                    }
                    for (int k = 0; k < 3; k += 1) { out_v(base_offset_v + num_vertices_per_geom - 1, k) = vi[k]; }

                    for (int j = 0; j < num_faces_per_geom; j += 1) {
                        const int findex0 = num_vertices_per_geom - 1;
                        const int findex1 = j;
                        const int findex2 = (j + 1) % (num_vertices_per_geom - 1);
                        out_f(base_offset_f + j, 0) = base_offset_v + findex0;
                        out_f(base_offset_f + j, 1) = base_offset_v + findex1;
                        out_f(base_offset_f + j, 2) = base_offset_v + findex2;
                    }
                    break;
                }
                case GeometryType::SPHERE: {
                    make_sphere_geometry(geometry_subdivisions_1, geometry_subdivisions_2,
                                         geometry_radius(i, 0),
                                         base_offset_f, base_offset_v,
                                         v(i, 0), v(i, 1), v(i, 2),
                                         out_v, out_f);
                    // throw std::runtime_error("Not implemented.");
                    break;
                }
            }
        }

        return num_faces_per_geom;
    }

    template <typename TRO, typename TRD>
    void trace_rays_point_cloud(const TRO& ray_o, const TRD& ray_d, double ray_near, double ray_far,
                                bool use_single_ray_origin,
                                int num_faces_per_geometry,
                                const EmbreeIntersector& isector,
                                Eigen::Matrix<typename TRO::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_t,
                                Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& out_pid) {
        for (int i = 0; i < ray_d.rows(); i++) {
            Eigen::RowVector3f o_i;
            if (use_single_ray_origin) {
                o_i = Eigen::RowVector3f((float) ray_o(0, 0), (float) ray_o(1, 0), (float) ray_o(2, 0));
            } else {
                o_i = Eigen::RowVector3f((float) ray_o(i, 0), (float) ray_o(i, 1), (float) ray_o(i, 2));
            }
            Eigen::RowVector3f d_i((float) ray_d(i, 0), (float) ray_d(i, 1), (float) ray_d(i, 2));
            Hit hit;

            bool is_hit = isector.intersectRay(o_i, d_i, hit, ray_near, ray_far);
            if (is_hit) {
                out_pid(i, 0) = (int) hit.id / num_faces_per_geometry;
                out_t(i, 0) = (typename TRO::Scalar) hit.t;
            } else {
                out_pid(i, 0) = -1;
                out_t(i, 0) = std::numeric_limits<typename TRO::Scalar>::infinity();
            }
        }
    }
}

npe_function(_populate_ray_point_cloud_intersector_internal)
npe_arg(v, dense_float, dense_double)
npe_arg(n, npe_matches(v))
npe_arg(geometry_type, std::string)
npe_arg(geometry_radius, npe_matches(v))
npe_default_arg(geometry_subdivisions_1, int, 4)
npe_default_arg(geometry_subdivisions_2, int, 4)
npe_arg(isector, std::shared_ptr<EmbreeIntersector>)
npe_begin_code()
    using MatrixI = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MatrixF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    GeometryType geom_type = validate_point_geometry(v, n, geometry_radius,
                                                     geometry_subdivisions_1,
                                                     geometry_subdivisions_2,
                                                     geometry_type);
    MatrixF geom_vertices;
    MatrixI geom_faces;

    int num_faces_per_geometry = generate_splat_geometry(geom_type, geometry_subdivisions_1, geometry_subdivisions_2,
                                                         v, n, geometry_radius, geom_vertices, geom_faces);
    isector->init(geom_vertices, geom_faces, true /*is_static*/);

    return num_faces_per_geometry;
npe_end_code()


npe_function(_intersect_ray_point_cloud_intersector_internal)
npe_arg(ray_o, dense_float, dense_double)
npe_arg(ray_d, npe_matches(ray_o))
npe_arg(isector, std::shared_ptr<EmbreeIntersector>)
npe_arg(num_faces_per_geometry, int)
npe_default_arg(ray_near, double, 0.0)
npe_default_arg(ray_far, double, std::numeric_limits<double>::infinity())
npe_begin_code()
{
    using MatrixI = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MatrixFOut = Eigen::Matrix<npe_Scalar_ray_o, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    bool use_single_ray_origin = validate_rays(ray_o, ray_d);

    MatrixI ret_pid(ray_d.rows(), 1);
    MatrixFOut ret_t(ray_d.rows(), 1);

    trace_rays_point_cloud(ray_o, ray_d, ray_near, ray_far,
                           use_single_ray_origin, num_faces_per_geometry,
                           *isector, ret_t, ret_pid);

    return std::make_tuple(npe::move(ret_pid), npe::move(ret_t));
}
npe_end_code()


npe_function(point_cloud_splatting_geometry_internal_)
npe_arg(v, dense_float, dense_double)
npe_arg(n, npe_matches(v))
npe_arg(geometry_type, std::string)
npe_arg(geometry_radius, npe_matches(v))
npe_default_arg(geometry_subdivisions_1, int, 4)
npe_default_arg(geometry_subdivisions_2, int, 4)
npe_begin_code()
{
    GeometryType geom_type = geometry_type_from_string(geometry_type);
    Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> geom_vertices;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> geom_faces;

    generate_splat_geometry(geom_type, geometry_subdivisions_1, geometry_subdivisions_2, v, n, geometry_radius,
                            geom_vertices, geom_faces);


    return std::make_tuple(npe::move(geom_vertices), npe::move(geom_faces));
}
npe_end_code()


npe_function(ray_point_cloud_intersection_internal_)
npe_arg(v, dense_float, dense_double)
npe_arg(n, npe_matches(v))
npe_arg(ray_o, npe_matches(v))
npe_arg(ray_d, npe_matches(v))
npe_arg(geometry_type, std::string)
npe_arg(geometry_radius, npe_matches(v))
npe_default_arg(geometry_subdivisions_1, int, 4)
npe_default_arg(geometry_subdivisions_2, int, 4)
npe_default_arg(ray_near, double, 0.0)
npe_default_arg(ray_far, double, std::numeric_limits<double>::infinity())
npe_begin_code()
    {
        using MatrixI = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using MatrixF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using MatrixFOut = Eigen::Matrix<npe_Scalar_v, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        bool use_single_ray_origin = validate_rays(ray_o, ray_d);
        GeometryType geom_type = validate_point_geometry(v, n, geometry_radius,
                                                         geometry_subdivisions_1,
                                                         geometry_subdivisions_2,
                                                         geometry_type);

        MatrixF geom_vertices;
        MatrixI geom_faces;

        int num_faces_per_geometry = generate_splat_geometry(geom_type, geometry_subdivisions_1, geometry_subdivisions_2,
                                                             v, n, geometry_radius, geom_vertices, geom_faces);

        EmbreeIntersector isector;
        isector.init(geom_vertices, geom_faces, true /*is_static*/);

        MatrixI ret_pid(ray_d.rows(), 1);
        MatrixFOut ret_t(ray_d.rows(), 1);

        trace_rays_point_cloud(ray_o, ray_d, ray_near, ray_far,
                               use_single_ray_origin, num_faces_per_geometry,
                               isector, ret_t, ret_pid);

        return std::make_tuple(npe::move(ret_pid), npe::move(ret_t));
}
npe_end_code()
