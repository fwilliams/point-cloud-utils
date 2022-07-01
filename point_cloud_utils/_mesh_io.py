import numpy as np
from contextlib import contextmanager
import os


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
            colors: [V, 4]-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] or [0, 255] (or None)
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
            self.custom_attributes = dict()
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

        def _set_color(self, key, clr_array):
            if not isinstance(self.colors, np.ndarray):
                self.colors = np.ndarray([clr_array.shape[0], 4], dtype=clr_array.dtype)
            else:
                assert clr_array.dtype == self.colors.dtype
            self.colors = self.colors.reshape([clr_array.shape[0], 4])
            if key == "alpha":
                self.colors[:, -1] = np.squeeze(clr_array)
            elif key == "colors":
                if clr_array.shape[1] == 3:
                    self.colors[:, :-1] = clr_array
                elif clr_array.shape[1] == 4:
                    self.colors = clr_array
                else:
                    raise ValueError("Invalid shape for vertex colors must be "
                                     "(n, 3) or (n, 4) but got " + str(clr_array.shape))
            else:
                assert False

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
            colors: [F, 4]-shaped numpy array of per-face RBGA colors in [0.0, 1.0] or [0, 255] (or None)
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
            self.custom_attributes = dict()
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

        def _set_color(self, key, clr_array):
            if not isinstance(self.colors, np.ndarray):
                self.colors = np.ndarray([clr_array.shape[0], 4], dtype=clr_array.dtype)
            else:
                assert clr_array.dtype == self.colors.dtype
            self.colors = self.colors.reshape([clr_array.shape[0], 4])
            if key == "alpha":
                self.colors[:, -1] = np.squeeze(clr_array)
            elif key == "colors":
                if clr_array.shape[1] == 3:
                    self.colors[:, :-1] = clr_array
                elif clr_array.shape[1] == 4:
                    self.colors = clr_array
                else:
                    raise ValueError("Invalid shape for face colors must be "
                                     "(n, 3) or (n, 4) but got " + str(clr_array.shape))
            else:
                assert False

        def _set_empty_to_none(self):
            for k, v in self.__dict__.items():
                if isinstance(v, np.ndarray):
                    if v.size <= 0:
                        setattr(self, k, None)

    def __init__(self, filename=None, dtype=np.float64):
        self.vertex_data = TriangleMesh.VertexData()
        self.face_data = TriangleMesh.FaceData()
        self.textures = []
        self.normal_maps = []
        if filename is not None:
            self.load(filename, dtype=dtype)

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

    def save(self, filename, dtype=np.float32):
        from ._pcu_internal import save_mesh_internal
        self.vertex_data._reset_if_none()
        self.face_data._reset_if_none()

        # Handle RBG colors by just concatenating alpha=1
        vcolors = self.vertex_data.colors
        if vcolors.shape[-1] == 3 and len(vcolors.shape) == 2:
            vcolors = np.concatenate([np.ascontiguousarray(self.vertex_data.colors),
                                      np.ones([vcolors.shape[0], 1], dtype=vcolors.dtype)], axis=-1)

        fcolors = self.face_data.colors
        if fcolors.shape[-1] == 3 and len(fcolors.shape) == 2:
            # vcglib face colors have an inverted alpha channel (so 1 = transparent, 0 = opaque)
            fcolors = np.concatenate([np.ascontiguousarray(fcolors),
                                      np.zeros([fcolors.shape[0], 1], dtype=self.face_data.colors.dtype)], axis=-1)
        elif fcolors.shape[-1] == 4 and len(fcolors.shape) == 2:
            fcolors[:, -1] = 1.0 - fcolors[:, -1]

        wcolors = self.face_data.wedge_colors
        if wcolors.shape[-1] == 3 and len(wcolors.shape) == 3:
            wcolors = np.concatenate([np.ascontiguousarray(wcolors),
                                      np.ones([wcolors.shape[0], wcolors.shape[1], 1], dtype=wcolors.dtype)], axis=-1)

        if fcolors.shape[0] > 0:
            if fcolors.dtype == np.uint8:
                fcolors = fcolors.astype(dtype) / 255.0
            if fcolors.max() > 1.0 or fcolors.min() < 0.0:
                raise ValueError("Invalid values for face colors, must be between 0 and 1 (inclusive)")
        if vcolors.shape[0] > 0:
            if vcolors.dtype == np.uint8:
                vcolors = vcolors.astype(dtype) / 255.0
            if vcolors.max() > 1.0 or vcolors.min() < 0.0:
                raise ValueError("Invalid values for vertex colors, must be between 0 and 1 (inclusive)")
        if wcolors.shape[0] > 0:
            if wcolors.dtype == np.uint8:
                wcolors = wcolors.astype(dtype) / 255.0
            if wcolors.max() > 1.0 or wcolors.min() < 0.0:
                raise ValueError("Invalid values for wedge colors, must be between 0 and 1 (inclusive)")

        save_mesh_internal(filename,
                           np.ascontiguousarray(self.vertex_data.positions.astype(dtype)),
                           np.ascontiguousarray(self.vertex_data.normals.astype(dtype)),
                           np.ascontiguousarray(self.vertex_data.texcoords.astype(dtype)),
                           np.ascontiguousarray(vcolors.astype(dtype)),
                           np.ascontiguousarray(self.vertex_data.quality.astype(dtype)),
                           np.ascontiguousarray(self.vertex_data.radius.astype(dtype)),
                           np.ascontiguousarray(self.vertex_data.tex_ids.astype(np.int32)),
                           np.ascontiguousarray(self.vertex_data.flags.astype(np.int32)),

                           np.ascontiguousarray(self.face_data.vertex_ids.astype(np.int32)),
                           np.ascontiguousarray(self.face_data.normals.astype(dtype)),
                           np.ascontiguousarray(fcolors.astype(dtype)),
                           np.ascontiguousarray(self.face_data.quality.astype(dtype)),
                           np.ascontiguousarray(self.face_data.flags.astype(np.int32)),

                           np.ascontiguousarray(wcolors.astype(dtype)),
                           np.ascontiguousarray(self.face_data.wedge_normals.astype(dtype)),
                           np.ascontiguousarray(self.face_data.wedge_texcoords.astype(dtype)),
                           np.ascontiguousarray(self.face_data.wedge_tex_ids.astype(np.int32)),
                           {k: np.ascontiguousarray(v) if isinstance(v, np.ndarray) else v for (k, v) in self.vertex_data.custom_attributes.items()},
                           {k: np.ascontiguousarray(v) if isinstance(v, np.ndarray) else v for (k, v) in self.face_data.custom_attributes.items()},
                           self.textures, self.normal_maps, dtype, np.int32)
        self.vertex_data._set_empty_to_none()
        self.face_data._set_empty_to_none()

    def load(self, filename, dtype=np.float64):
        from ._pcu_internal import load_mesh_internal
        if not os.path.exists(filename):
            raise FileNotFoundError("Invalid path " + filename + " does not exist")
        mesh_root_path = os.path.dirname(filename)
        if mesh_root_path.strip() == '':
            mesh_root_path = '.'
        mesh_filename = os.path.basename(filename)

        @contextmanager
        def pushd(new_dir):
            previous_dir = os.getcwd()
            os.chdir(new_dir)
            try:
                yield
            finally:
                os.chdir(previous_dir)

        with pushd(mesh_root_path):
            mesh_dict = load_mesh_internal(mesh_filename, dtype)

        self.textures = mesh_dict["textures"]
        self.normal_maps = mesh_dict["normal_maps"]
        vret = mesh_dict["vertex_data"]
        fret = mesh_dict["face_data"]
        for k, v in vret.items():
            if k in ("alpha", "colors"):
                self.vertex_data._set_color(k, v)
            else:
                assert hasattr(self.vertex_data, k), "vertex_data doesn't have attribute " + str(k)
                setattr(self.vertex_data, k, v)
        for k, v in fret.items():
            if k in ("alpha", "colors"):
                self.face_data._set_color(k, v)
            else:
                assert hasattr(self.face_data, k), "face_data doesn't have attribute " + str(k)
                setattr(self.face_data, k, v)
        self.vertex_data._set_empty_to_none()
        self.face_data._set_empty_to_none()


def save_triangle_mesh(filename, v, f=None,
                       vn=None, vt=None, vc=None, vq=None, vr=None, vti=None, vflags=None,
                       fn=None, fc=None, fq=None, fflags=None, wc=None, wn=None, wt=None, wti=None,
                       textures=[], normal_maps=[], dtype=np.float32):
    """
    Save a triangle mesh to a file with various per-vertex, per-face, and per-wedge attributes. Each argument (except v)
    is optional and can be None.

    Parameters
    ----------
    filename    : Path to the mesh to save. The type of file will be determined from the file extension.
    v           : [V, 3]-shaped numpy array of per-vertex positions
    f           : [F, 3]-shaped numpy array of integer face indices into TrianglMesh.vertex_data.positions (or None)
    vn          : [V, 3]-shaped numpy array of per-vertex normals (or None)
    vt          : [V, 2]-shaped numpy array of per-vertex uv coordinates (or None)
    vc          : [V, 4]-shaped numpy array of per-vertex RBGA colors in [0.0, 1.0] (or None)
    vq          : [V,]-shaped numpy array of per-vertex quality measures (or None)
    vr          : [V,]-shaped numpy array of per-vertex curvature radii (or None)
    vti         : [V,]-shaped numpy array of integer indices into TriangleMesh.textures indicating which texture to
                  use at this vertex (or None)
    vflags      : [V,]-shaped numpy array of 32-bit integer flags per vertex (or None)
    fn          : [F, 3]-shaped numpy array of per-face normals (or None)
    fc          : [F, 4]-shaped numpy array of per-face RBGA colors in [0.0, 1.0] (or None)
    fq          : [F,]-shaped numpy array of per-face quality measures (or None)
    fflags      : [F,]-shaped numpy array of 32-bit integer flags per face (or None)
    wc          : [F, 3, 4]-shaped numpy array of per-wedge RBGA colors in [0.0, 1.0] (or None)
    wn          : [F, 3, 3]-shaped numpy array of per-wedge normals (or None)
    wt          : [F, 3, 2]-shaped numpy array of per-wedge] uv coordinates (or None)
    wti         : [F, 3]-shaped numpy array of integer indices into TriangleMesh.textures indicating which
    textures    : A list of paths to texture image files for this mesh
    normal_maps : A list of paths to texture image files for this mesh

    Returns
    -------
    None
    """
    mesh = TriangleMesh()
    mesh.vertex_data.positions = v
    mesh.vertex_data.normals = vn
    mesh.vertex_data.texcoords = vt
    mesh.vertex_data.colors = vc
    mesh.vertex_data.quality = vq
    mesh.vertex_data.radius = vr
    mesh.vertex_data.tex_ids = vti
    mesh.vertex_data.flags = vflags

    mesh.face_data.vertex_ids = f
    mesh.face_data.normals = fn
    mesh.face_data.colors = fc
    mesh.face_data.quality = fq
    mesh.face_data.flags = fflags

    mesh.face_data.wedge_colors = wc
    mesh.face_data.wedge_normals = wn
    mesh.face_data.wedge_texcoords = wt
    mesh.face_data.wedge_tex_ids = wti

    mesh.textures = textures
    mesh.normal_maps = normal_maps

    mesh.save(filename, dtype=dtype)


def save_mesh_v(filename, v, dtype=np.float32):
    save_triangle_mesh(filename, v=v, dtype=dtype)


def save_mesh_vf(filename, v, f, dtype=np.float32):
    save_triangle_mesh(filename, v=v, f=f, dtype=dtype)


def save_mesh_vn(filename, v, n, dtype=np.float32):
    save_triangle_mesh(filename, v=v, vn=n, dtype=dtype)


def save_mesh_vc(filename, v, c, dtype=np.float32):
    save_triangle_mesh(filename, v=v, vc=c, dtype=dtype)


def save_mesh_vnc(filename, v, n, c, dtype=np.float32):
    save_triangle_mesh(filename, v=v, vn=n, vc=c, dtype=dtype)


def save_mesh_vfn(filename, v, f, n, dtype=np.float32):
    save_triangle_mesh(filename, v=v, f=f, vn=n, dtype=dtype)


def save_mesh_vfnc(filename, v, f, n, c, dtype=np.float32):
    save_triangle_mesh(filename, v=v, f=f, vn=n, vc=c, dtype=dtype)


def load_triangle_mesh(filename, dtype=np.float64):
    ret = TriangleMesh()
    ret.load(filename, dtype=dtype)
    return ret


def load_mesh_v(filename, dtype=float):
    return load_triangle_mesh(filename, dtype=dtype).v


def load_mesh_vf(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.f


def load_mesh_vn(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.vn


def load_mesh_vc(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.vc


def load_mesh_vnc(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.vn, ret.vc


def load_mesh_vfn(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.f, ret.vn


def load_mesh_vfnc(filename, dtype=float):
    ret = load_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.f, ret.vn, ret.vc

