import os
from random import choice
from compas.geometry import Polyline, Polygon
from compas.datastructures import Mesh, mesh_thicken
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput

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
plate_thickness = (10*units.cm).to_base_units().magnitude
plate = mesh_thicken(mesh, plate_thickness)

# ==============================================================================
# Select random internal vertices for the point, line and area load application
# ==============================================================================
vertex = choice(list(set(mesh.vertices()) - set(mesh.vertices_on_boundary())))
poa = mesh.vertex_coordinates(vertex)

points = [choice(list(set(mesh.vertices()) - set(mesh.vertices_on_boundary()))) for _ in range(3)]
polyline = Polyline([mesh.vertex_coordinates(p) for p in points])

face = choice(list(set(mesh.faces())))
polygon = Polygon([[p[0], p[1], p[2]+plate_thickness/2] for p in mesh.face_coordinates(face)])

# ==============================================================================
# GMSH model
# ==============================================================================
targetlength = 400
model = MeshModel.from_mesh(plate, targetlength=targetlength)
model.heal()
model.refine_mesh()
model.generate_mesh(3)
solid_mesh = model.mesh_to_compas()

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='mesh_refine')
# Define some properties
mat = ElasticIsotropic(E=(210*units.GPa).to_base_units().magnitude,
                       v=0.2,
                       density=(7800*units("kg/m**3")).to_base_units().magnitude)
sec = SolidSection(material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec)
prt.ndf = 3  # this is needed for the opensees FourNodeTetrahedron model
prt._discretized_boundary_mesh = solid_mesh
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in mesh.vertices_where({'vertex_degree': 2}):
    location = mesh.vertex_coordinates(vertex)
    mdl.add_fix_bc(nodes=prt.find_nodes_by_location(location, distance=150))

# mdl.summary()
# mdl.show()

# Set-up the problem
prb = Problem('01_mesh_refine')
mdl.add_problem(problem=prb)

# Initialize a step
stp = StaticStep()
prb.add_step(stp)

# Add the loads
stp.add_point_load(points=[poa], z=-(300*units.N).to_base_units().magnitude)
stp.add_line_load(polyline, z=(-1*units.kN).to_base_units().magnitude)
stp.add_planar_area_load(polygon, z=(-3*units.kN).to_base_units().magnitude)

# Ask for field outputs
stp.add_output(FieldOutput(node_outputs=['U']))

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Show Results
prb.show_displacements(draw_loads=5, draw_bcs=500)
