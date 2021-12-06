import numpy as np

class TriangleMesh:
    """
    A lightweight container class representing a triangle mesh with attributes stored at each vertex, wedge, and face,
    where:
        - A vertex is a 3D position
        - A wedge is a vertex and its two adjacent edges
        - A face is a triangle connecting 3 vertices (Each triangle has 3 vertices and 3 wedges).

    The data in this class is represented as follows:

    TriangleMesh:
        vertex_data:
            positions: [V, 3]-shaped numpy array of per-vertex positions
            normals: [V, 3]-shaped numpy array of per-vertex normals (or None)
            texcoords: [V, 2]-shaped numpy array of per-vertex uv coordinates (or None)
            tex_ids: [V,]-shaped numpy array of integer indices into TriangleMesh.textures indicating which texture to
                     use at this vertex (or None)
            colors: [V, 4]-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] (or None)
            radius: [V,]-shaped numpy array of per-vertex curvature radii (or None)
            quality: [V,]-shaped numpy array of per-vertex quality measures (or None)
            flags: [V,]-shaped numpy array of 32-bit integer flags per vertex (or None)
        face_data:
            vertex_ids: [F, 3]-shaped numpy array of integer face indices into TrianglMesh.vertex_data.positions
            normals: [F, 3]-shaped numpy array of per-face normals (or None)
            colors: [F, 4]-shaped numpy array of per-face RBGA colors in [0.0, 1.0] (or None)
            quality: [F,]-shaped numpy array of per-face quality measures (or None)
            flags: [F,]-shaped numpy array of 32-bit integer flags per face (or None)

            wedge_colors: [F, 3, 4]-shaped numpy array of per-wedge RBGA colors in [0.0, 1.0] (or None)
            wedge_normals: [F, 3, 3]-shaped numpy array of per-wedge normals (or None)
            wedge_texcoords: [F, 3, 2]-shaped numpy array of per-wedge] uv coordinates (or None)
            wedge_tex_ids: [F, 3]-shaped numpy array of integer indices into TriangleMesh.textures indicating which
                           texture to use at this wedge (or None)
        textures: A list of paths to texture image files for this mesh
        normal_maps: A list of paths to texture image files for this mesh
    """
    class VertexData:
        """
        A lightweight container class representing per-vertex information within a TriangleMesh.

        The data in this class is represented as follows:

        VertexData:
            positions: [V, 3]-shaped numpy array of per-vertex positions
            normals: [V, 3]-shaped numpy array of per-vertex normals (or None)
            texcoords: [V, 2]-shaped numpy array of per-vertex uv coordinates (or None)
            tex_ids: [V,]-shaped numpy array of integer indices into TriangleMesh.textures indicating which texture to
                     use at this vertex (or None)
            colors: [V, 4]-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] (or None)
            radius: [V,]-shaped numpy array of per-vertex curvature radii (or None)
            quality: [V,]-shaped numpy array of per-vertex quality measures (or None)
            flags: [V,]-shaped numpy array of 32-bit integer flags per vertex (or None)
        """

        def __init__(self):
            self.positions = np.zeros([0, 3])
            self.normals = np.zeros([0, 3])
            self.texcoords = np.zeros([0, 2])
            self.colors = np.zeros([0, 4])
            self.quality = np.zeros([0])
            self.radius = np.zeros([0])
            self.tex_ids = np.zeros([0], dtype=int)
            self.flags = np.zeros([0], dtype=int)
            self._set_empty_to_none()

        def _reset_if_none(self):
            if self.positions is None:
                self.positions = np.zeros([0, 3])
            if self.normals is None:
                self.normals = np.zeros([0, 3])
            if self.texcoords is None:
                self.texcoords = np.zeros([0, 2])
            if self.colors is None:
                self.colors = np.zeros([0, 4])
            if self.quality is None:
                self.quality = np.zeros([0])
            if self.radius is None:
                self.radius = np.zeros([0])
            if self.tex_ids is None:
                self.tex_ids = np.zeros([0], dtype=int)
            if self.flags is None:
                self.flags = np.zeros([0], dtype=int)

        def _set_empty_to_none(self):
            for k, v in self.__dict__.items():
                if isinstance(v, np.ndarray):
                    if v.size <= 0:
                        setattr(self, k, None)

    class FaceData:
        """
        A lightweight container class representing per-face information within a TriangleMesh.

        The data in this class is represented as follows:

        FaceData:
            vertex_ids: [F, 3]-shaped numpy array of integer face indices into TrianglMesh.vertex_data.positions
            normals: [F, 3]-shaped numpy array of per-face normals (or None)
            colors: [F, 4]-shaped numpy array of per-face RBGA colors in [0.0, 1.0] (or None)
            quality: [F,]-shaped numpy array of per-face quality measures (or None)
            flags: [F,]-shaped numpy array of 32-bit integer flags per face (or None)

            wedge_colors: [F, 3, 4]-shaped numpy array of per-wedge RBGA colors in [0.0, 1.0] (or None)
            wedge_normals: [F, 3, 3]-shaped numpy array of per-wedge normals (or None)
            wedge_texcoords: [F, 3, 2]-shaped numpy array of per-wedge] uv coordinates (or None)
            wedge_tex_ids: [F, 3]-shaped numpy array of integer indices into TriangleMesh.textures indicating which
                           texture to use at this wedge (or None)
        """

        def __init__(self):
            self.vertex_ids = np.zeros([0, 3], dtype=int)
            self.normals = np.zeros([0, 3])
            self.colors = np.zeros([0, 4])
            self.quality = np.zeros([0])
            self.flags = np.zeros([0], dtype=int)

            self.wedge_colors = np.zeros([0, 3, 4])
            self.wedge_normals = np.zeros([0, 3, 3])
            self.wedge_texcoords = np.zeros([0, 3, 2])
            self.wedge_tex_ids = np.zeros([0, 3], dtype=int)
            self._set_empty_to_none()

        def _reset_if_none(self):
            if self.vertex_ids is None:
                self.vertex_ids = np.zeros([0, 3], dtype=int)
            if self.normals is None:
                self.normals = np.zeros([0, 3])
            if self.colors is None:
                self.colors = np.zeros([0, 4])
            if self.quality is None:
                self.quality = np.zeros([0])
            if self.flags is None:
                self.flags = np.zeros([0], dtype=int)

            if self.wedge_colors is None:
                self.wedge_colors = np.zeros([0, 3, 4])
            if self.wedge_normals is None:
                self.wedge_normals = np.zeros([0, 3, 3])
            if self.wedge_texcoords is None:
                self.wedge_texcoords = np.zeros([0, 3, 2])
            if self.wedge_tex_ids is None:
                self.wedge_tex_ids = np.zeros([0, 3], dtype=int)

        def _set_empty_to_none(self):
            for k, v in self.__dict__.items():
                if isinstance(v, np.ndarray):
                    if v.size <= 0:
                        setattr(self, k, None)

    def __init__(self, dtype=np.float64):
        self.vertex_data = TriangleMesh.VertexData()
        self.face_data = TriangleMesh.FaceData()
        self.textures = []
        self.normal_maps = []

    @property
    def v(self):
        return self.vertex_data.positions

    @property
    def f(self):
        return self.face_data.vertex_ids

    @property
    def vn(self):
        return self.vertex_data.normals

    @property
    def vc(self):
        return self.vertex_data.colors

    @property
    def fn(self):
        return self.vertex_data.normals

    @property
    def fc(self):
        return self.face_data.colors
