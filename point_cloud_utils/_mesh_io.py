import numpy as np


class TriangleMesh:
    class VertexData:
        def __init__(self):
            self.positions = np.zeros([0, 3])
            self.normals = np.zeros([0, 3])
            self.texcoords = np.zeros([0, 2])
            self.colors = np.zeros([0, 4])
            self.quality = np.zeros([0])
            self.radius = np.zeros([0])
            self.tex_ids = np.zeros([0], dtype=int)
            self.flags = np.zeros([0], dtype=int)

    class FaceData:
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

    def __init__(self):
        self.vertex_data = TriangleMesh.VertexData()
        self.face_data = TriangleMesh.FaceData()
        self.textures = []
        self.normal_maps = []

        # Shortcuts to commonly accessed attributes
        self.v = self.vertex_data.positions
        self.f = self.face_data.vertex_ids
        self.vn = self.vertex_data.normals
        self.fn = self.face_data.normals


def read_triangle_mesh(filename, dtype=float):
    from ._pcu_internal import load_mesh
    mesh_dict = load_mesh(filename, dtype)
    ret = TriangleMesh()
    ret.textures = mesh_dict["textures"]
    ret.normal_maps = mesh_dict["normal_maps"]
    vret = mesh_dict["vertex_data"]
    fret = mesh_dict["face_data"]
    for k, v in vret.items():
        assert hasattr(ret.vertex_data, k)
        setattr(ret.vertex_data, k, v)
    for k, v in fret.items():
        assert hasattr(ret.face_data, k)
        setattr(ret.face_data, k, v)
    ret.v = ret.vertex_data.positions
    ret.f = ret.face_data.vertex_ids
    ret.vn = ret.vertex_data.normals
    ret.fn = ret.face_data.normals

    return ret


def read_mesh_v(filename, dtype=float):
    return read_triangle_mesh(filename, dtype=dtype).v


def read_mesh_vf(filename, dtype=float):
    ret = read_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.f


def read_mesh_vfn(filename, dtype=float):
    ret = read_triangle_mesh(filename, dtype=dtype)
    return ret.v, ret.f, ret.vn