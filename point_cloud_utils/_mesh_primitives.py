import numpy as np

# These are all writen by ChatGPT lol :)

def cylinder_mesh(radius, height, sides=16, top=True, bottom=True):
    """
    Generate a triangle mesh for a cylinder centered at the origin aligned with the y-axis.
    
    Args:
        radius (float) : Radius of the cylinder.
        height (float) : Height of the cylinder.
        sides (int) : Number of sides (segments) for the cylinder (default = 16).
        top (bool) : Include top face if True (default = True).
        bottom (bool) : Include bottom face if True (default = True).
        
    Returns:
        vertices (np.ndarray) : List of vertex coordinates (shape = (v, 3)).
        faces (np.ndarray) : List of faces (triplets of vertex indices with shape = (f, 3)).
    """
    
    if not isinstance(radius, (float, int)):
        raise ValueError("radius must be a number")
    if not isinstance(height, (float, int)):
        raise ValueError("height must be a number")
    if not isinstance(sides, int):
        raise ValueError("sides must be an integer")
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if height <= 0.0:
        raise ValueError("height must be positive")
    # Calculate the angle between each side of the cylinder
    angle_step = 2 * np.pi / sides
    
    # Generate vertices for the sides of the cylinder
    vertices = []
    for i in range(sides):
        x = radius * np.cos(i * angle_step)
        y = -height / 2
        z = radius * np.sin(i * angle_step)
        vertices.append((x, y, z))
        vertices.append((x, y + height, z))
    
    # Generate vertices for the top and bottom if required
    if top:
        vertices.append((0, -height / 2, 0))
    if bottom:
        vertices.append((0, height / 2, 0))
    
    # Generate faces for the sides of the cylinder
    faces = []
    for i in range(sides):
        v1 = 2 * i
        v2 = (2 * i + 2) % (2 * sides)
        v3 = (2 * i + 1) % (2 * sides)
        v4 = (2 * i + 3) % (2 * sides)
        faces.append((v1, v3, v2))
        faces.append((v4, v2, v3))
    
    # Generate faces for the top and bottom if required
    if top:
        top_vertex = len(vertices) - 2
        for i in range(sides):
            v1 = top_vertex
            v2 = (2 * i) % (2 * sides)
            v3 = (2 * (i + 1)) % (2 * sides)
            faces.append((v1, v2, v3))
    
    if bottom:
        bottom_vertex = len(vertices) - 1
        for i in range(sides):
            v1 = bottom_vertex
            v2 = (2 * i + 1) % (2 * sides)
            v3 = (2 * (i + 1) + 1) % (2 * sides)
            faces.append((v1, v3, v2))
    
    return np.array(vertices), np.array(faces)



def cube_mesh():
    """
    Generate a triangle mesh for a unit cube centered at the origin.

    Returns:
        vertices (np.ndarray) : Mesh vertices as an (n, 3)-shaped NumPy array
        faces (np.ndarray) : Mesh face indices as an (f, 3)-shaped NumPy array
    """
    size = 0.5
    vertices = [
        [-size, -size, -size],
        [-size, -size, size],
        [-size, size, -size],
        [-size, size, size],
        [size, -size, -size],
        [size, -size, size],
        [size, size, -size],
        [size, size, size]
    ]
    
    # Define the indices that form the triangles of the cube's faces
    indices = [
        [0, 1, 2], [1, 3, 2], # Front face
        [4, 6, 5], [5, 6, 7], # Back face
        [0, 2, 4], [4, 2, 6], # Left face
        [1, 5, 3], [5, 7, 3], # Right face
        [2, 3, 6], [3, 7, 6], # Top face
        [0, 4, 1], [4, 5, 1]  # Bottom face
    ]
    
    vertices = np.array(vertices, dtype=np.float32)
    indices = np.array(indices, dtype=np.uint32)
    
    return vertices, indices


def sphere_mesh(subdivisions=3):
    """
    Generate a triangle mesh approximating a unit sphere centered at the origin by subdividing an isocahedron.
    
    Args:
        subdivisions (int) : Number of times to subdivide an isocahedron to get the mesh (default = 3).

    Returns:
        vertices (np.ndarray) : Mesh vertices as an (n, 3)-shaped NumPy array
        faces (np.ndarray) : Mesh face indices as an (f, 3)-shaped NumPy array
    """

    if not isinstance(subdivisions, int):
        raise ValueError("Invalid type for subdivisions, must be an integer")

    # Define the initial icosahedron vertices
    t = (1.0 + np.sqrt(5.0)) / 2.0
    
    vertices = np.array([
        [-1, t, 0],
        [1, t, 0],
        [-1, -t, 0],
        [1, -t, 0],
        [0, -1, t],
        [0, 1, t],
        [0, -1, -t],
        [0, 1, -t],
        [t, 0, -1],
        [t, 0, 1],
        [-t, 0, -1],
        [-t, 0, 1]
    ], dtype=np.float32)

    # Define the initial icosahedron triangles
    triangles = np.array([
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1]
    ], dtype=np.uint32)

    for _ in range(subdivisions):
        new_vertices = []
        new_triangles = []
        
        for triangle in triangles:
            v1 = vertices[triangle[0]]
            v2 = vertices[triangle[1]]
            v3 = vertices[triangle[2]]
            
            mid1 = (v1 + v2) / 2
            mid2 = (v2 + v3) / 2
            mid3 = (v3 + v1) / 2
            
            new_vertices.extend([v1, v2, v3, mid1, mid2, mid3])
            
            new_v1_idx = len(new_vertices) - 6
            new_v2_idx = len(new_vertices) - 5
            new_v3_idx = len(new_vertices) - 4
            new_mid1_idx = len(new_vertices) - 3
            new_mid2_idx = len(new_vertices) - 2
            new_mid3_idx = len(new_vertices) - 1
            
            new_triangles.extend([
                [new_v1_idx, new_mid1_idx, new_mid3_idx],
                [new_mid1_idx, new_v2_idx, new_mid2_idx],
                [new_mid3_idx, new_mid2_idx, new_v3_idx],
                [new_mid1_idx, new_mid2_idx, new_mid3_idx]
            ])
        
        vertices = np.array(new_vertices, dtype=np.float32)
        triangles = np.array(new_triangles, dtype=np.uint32)
    
    # Normalize the vertices to get the unit sphere
    vertices /= np.linalg.norm(vertices, axis=1, keepdims=True)
    
    return vertices, triangles


def cone_mesh(radius, height, sides=16, bottom=True):
    """
    Generate a triangle mesh for a cone centered at the origin and pointing up the y-axis.

    Args:
        radius (float): Radius of the cone's base.
        height (float): Height of the cone.
        sides (int): Number of segments used to approximate the circular base (default = 16).
        bottom (bool): Whether to include faces for the bottom/base of the code (default = True).

    Returns:
        vertices (np.ndarray) : Mesh vertices as an (n, 3)-shaped NumPy array
        faces (np.ndarray) : Mesh face indices as an (f, 3)-shaped NumPy array
    """
    vertices = []
    faces = []

    # Generate vertices for the cone sides
    for i in range(sides):
        angle = 2 * np.pi * i / sides
        x = radius * np.cos(angle)
        y = -height/2
        z = radius * np.sin(angle)
        vertices.append((x, y, z))

    # Apex of the cone
    vertices.append((0.0, height/2, 0.0))

    # Generate faces for the cone sides
    for i in range(sides):
        if i == sides - 1:
            faces.append((i, sides, 0))
        else:
            faces.append((i, sides, i + 1))

    # Generate faces for the base
    if bottom:
        base_center = len(vertices)
        vertices.append((0.0, 0.0, 0.0))
        for i in range(sides):
            if i == sides - 1:
                faces.append((i, 0, base_center))
            else:
                faces.append((i, i + 1, base_center))

    return np.array(vertices), np.array(faces)