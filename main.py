import numpy as np
from scipy.spatial import ConvexHull

# Change these three factors whatever you need
frequency_num = 3
decimal_accuracy = 5
radius = 7500

def generate_icosahedron():
    # Generates the 12 vertices of a unit icosahedron
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    verts = np.array([
        [-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
        [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
        [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1]
    ])
    verts /= np.linalg.norm(verts[0])  # Normalize to unit sphere
    return verts


def get_single_triangle(verts):
    # Gets one triangle from the icosahedron.
    hull = ConvexHull(verts)
    return hull.simplices[0]  # Return first face


def subdivide_triangle(verts, triangle, frequency=frequency_num):
    # Subdivides one triangle, tracks edges, and returns new vertices + edges.
    v1, v2, v3 = verts[triangle]  # Get original triangle vertices
    vertex_map = {}  # Maps (i, j) indices to vertex positions
    edge_set = set()  # Stores edges (as vertex index pairs)

    new_vertices = []

    # Create grid of subdivided points
    for i in range(frequency + 1):
        for j in range(frequency + 1 - i):
            t1, t2 = i / frequency, j / frequency
            new_point = (1 - t1 - t2) * v1 + t1 * v2 + t2 * v3  # Calculate new vertice by interpolation
            new_point /= np.linalg.norm(new_point)  # Project onto sphere
            vertex_map[(i, j)] = len(new_vertices)  # Store index
            new_vertices.append(tuple(new_point))

    # Define edges in subdivided structure
    for i in range(frequency + 1):
        for j in range(frequency + 1 - i):
            idx = vertex_map[(i, j)]

            # Horizontal edge (i+1, j)
            if (i + 1, j) in vertex_map:
                edge_set.add((idx, vertex_map[(i + 1, j)]))

            # Diagonal edge (i, j+1)
            if (i, j + 1) in vertex_map:
                edge_set.add((idx, vertex_map[(i, j + 1)]))

            # Triangle diagonal (i+1, j-1)
            if (i + 1, j - 1) in vertex_map:
                edge_set.add((idx, vertex_map[(i + 1, j - 1)]))

    return np.array(new_vertices), edge_set


def compute_chord_factors(vertices, edges):
    # Computes unique chord factors based on recorded edges of single face
    lengths = []
    for v1, v2 in edges:
        length = np.linalg.norm(np.array(vertices[v1]) - np.array(vertices[v2]))
        lengths.append(length)

    return sorted(set(np.round(lengths, decimal_accuracy)))  # Unique chord factors

def get_strut_lengths(chord_factors):
    # Compute strut lengths from chord factors
    strut_lengths = []

    for value in chord_factors:
        strut_lengths.append(value * 7500)

    return strut_lengths


# Generate and process the icosahedron
vertices = generate_icosahedron()
triangle = get_single_triangle(vertices)
subdivided_vertices, edges = subdivide_triangle(vertices, triangle)
chord_factors = compute_chord_factors(subdivided_vertices, edges)
strut_lengths = get_strut_lengths(chord_factors)


print("Chord Factors:", chord_factors)
print("Strut Lengths:", strut_lengths)

