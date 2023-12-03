from compas_view2.app import App
from compas_view2.shapes import Arrow

import compas
from compas.datastructures import Mesh
from compas.geometry import Point
import random

# # Create an initial planar mesh (as shown in the previous example)
# vertices = [
#     [0, 0, 0],
#     [10000, 0, 0],
#     [10000, 10000, 0],
#     [0, 10000, 0]
# ]
# faces = [
#     [0, 1, 2, 3]
# ]
# mesh = Mesh.from_vertices_and_faces(vertices, faces)


lx = 10000
ly = 10000
nx = 5
ny = 5
mesh = Mesh.from_meshgrid(lx, nx, ly, ny)
vertices = [mesh.vertex_coordinates(v) for v in mesh.vertices()]


# Define the number of random points to add
num_points = 1

# Define the bounding box of the mesh
bbox = mesh.bounding_box()


# Create a list of faces to delete
faces_to_delete = list(mesh.faces())

# Create a list of vertices to keep
vertices_to_keep = list(mesh.vertices())

# Remove the original faces
for face_key in faces_to_delete:
    mesh.delete_face(face_key)

# Add random points within the bounding box
for _ in range(num_points):
    x = random.uniform(bbox[0][0], bbox[1][0])
    y = random.uniform(bbox[0][1], bbox[1][1])
    z = 0  # For a planar mesh, set the z-coordinate to zero
    point = Point(x, y, z)
    mesh.add_vertex(x=x, y=y, z=z)  # No need to specify a key

# Create new faces connecting the random point to the existing vertices
random_point_index = len(vertices)  # Index of the random point
for i in range(len(vertices)):
    mesh.add_face([i, (i + 1) % len(vertices), random_point_index])

# Print the updated number of vertices and faces
print(f"Vertices (including random points): {len(list(mesh.vertices()))}")
print(f"Faces (including non-overlapping triangles): {len(list(mesh.faces()))}")



# ==============================================================================
# Viz
# ==============================================================================

viewer = App(width=1600, height=900)

viewer.view.camera.rz = 0
viewer.view.camera.rx = -55
viewer.view.camera.tx = -5
viewer.view.camera.ty = -2
viewer.view.camera.distance = 10

viewer.view.camera.target = [5000, 5000, 100]
viewer.view.camera.position = [5000, 0, 5000]
viewer.view.camera.near = 1
viewer.view.camera.far = 100000
viewer.view.camera.scale = 1000
viewer.view.grid.cell_size = 1000

viewer.add(mesh)


viewer.run()