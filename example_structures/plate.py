import os

from random import choice
from compas.datastructures import Mesh, mesh_thicken
from compas.utilities import geometric_key_xy
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
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
nx = 2
ny = 2
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

poa = choice(list(set(plate.vertices()) - set(plate.vertices_on_boundary())))
poa_coordinates = plate.vertex_coordinates(poa)

# ==============================================================================
# GMSH model
# ==============================================================================

model = MeshModel.from_mesh(plate, targetlength=1000)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# model.optimize_mesh(niter=100)
# model.recombine_mesh()

# ==============================================================================
# COMPAS mesh
# ==============================================================================

compas_mesh = model.mesh_to_compas()
lengths = [compas_mesh.edge_length(*edge) for edge in compas_mesh.edges()]
print('Min length: ', min(lengths))
print('Max length: ', max(lengths))

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=(210*units.GPa).to_base_units().magnitude, 
                       v=0.2, 
                       density=(7800*units("kg/m**3")).to_base_units().magnitude)
sec = ShellSection(material=mat, t=100)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec, implementation='shelldkgt')
# prt.ndf=6 # this is needed for the opensees FourNodeTetrahedron model
prt._discretized_boundary_mesh = compas_mesh
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in plate.vertices_where({'vertex_degree': 2}):
    location = plate.vertex_coordinates(vertex)
    mdl.add_fix_bc(nodes=prt.find_nodes_by_location(location, distance=150))

# mdl.summary()
# mdl.show()

# Initialize a step
stp = StaticStep()

# Add the load
pt = prt.find_closest_nodes_to_point(poa_coordinates, distance=150)
stp.add_point_load(nodes=pt,
                      z=-(10*units.kN).to_base_units().magnitude)

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
# mdl.analyse(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Show Results
# prb.show_displacements(show_deformed=100)
# prb.show_deformed()
