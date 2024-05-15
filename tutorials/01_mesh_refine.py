import os
from pathlib import Path

from random import choice
from compas.datastructures import Mesh
from compas.geometry import Point, Vector
from compas.utilities import geometric_key_xy
from compas_gmsh.models import MeshModel
from compas_view2.app import App
from compas_view2.shapes import Arrow

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, SolidSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (10*units.m).to_base_units().magnitude
ly = (10*units.m).to_base_units().magnitude
nx = 5
ny = 5
mesh = Mesh.from_meshgrid(lx, nx, ly, ny)
plate = mesh.thickened(100)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

poa = choice(list(set(mesh.vertices()) - set(mesh.vertices_on_boundary())))
poa_coordinates = mesh.vertex_coordinates(poa)

# ==============================================================================
# GMSH model
# ==============================================================================

model = MeshModel.from_mesh(plate, targetlength=500)

# refine mesh at the point of application of the load
model.mesh_targetlength_at_vertex(poa, 100)

# # refine mesh at the supports
# for vertex in mesh.vertices_where({'vertex_degree': 2}):
#     a = geometric_key_xy(mesh.vertex_coordinates(vertex))
#     for vertex in plate.vertices():
#         b = geometric_key_xy(plate.vertex_coordinates(vertex))
#         if a == b:
#             model.mesh_targetlength_at_vertex(vertex, 10)

model.heal()
model.refine_mesh()
model.generate_mesh(3)
# model.optimize_mesh(niter=100)
# model.recombine_mesh()

# ==============================================================================
# COMPAS mesh
# ==============================================================================

solid_mesh = model.mesh_to_compas()
lengths = [solid_mesh.edge_length(edge) for edge in solid_mesh.edges()]
print('Min length: ', min(lengths))
print('Max length: ', max(lengths))

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

viewer.add(solid_mesh)

poa = Point(* plate.vertex_coordinates(poa))
start = poa + Vector(0, 0, 1000)
vector = Vector(0, 0, -1000)
load = Arrow(start, vector, body_width=0.02)

viewer.add(poa, pointsize=5)
viewer.add(load, facecolor=(1, 0, 0))

viewer.run()

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='mesh_refine')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = SolidSection(material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec)
prt._boundary_mesh = plate  
prt.ndf=3 # this is needed for the opensees FourNodeTetrahedron model
prt._discretized_boundary_mesh = solid_mesh
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in mesh.vertices_where({'vertex_degree': 2}):
    location = mesh.vertex_coordinates(vertex)
    mdl.add_pin_bc(nodes=prt.find_nodes_around_point(location, distance=150))

# mdl.summary()

# Initialize a step
stp = StaticStep()

# Add the load
pt = prt.find_closest_nodes_to_point(poa_coordinates, distance=150)
stp.add_node_pattern(nodes=pt,
                      z=-10*units.kN,
                      load_case="LL")
stp.combination = LoadCombination.ULS()
# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S', 'SF'])
stp.add_output(fout)

# Set-up the problem
prb = Problem('01_mesh_refine')
prb.add_step(stp)
# prb.summary()
mdl.add_problem(problem=prb)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Show Results
# prb.show_displacements_contour(draw_loads=0.1, draw_bcs=400)
prb.show_deformed(scale_factor=5000, draw_loads=0.1, draw_bcs=0.3)