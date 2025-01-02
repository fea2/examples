import os

from random import choice
from pprint import pprint

from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel
from compas.colors import ColorMap, Color
from compas.geometry import Point, Line

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (3*units.m).to_base_units().magnitude
ly = (3*units.m).to_base_units().magnitude
nx = 10
ny = 10
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

poa = choice(list(set(plate.vertices()) - set(plate.vertices_on_boundary())))
poa_coordinates = plate.vertex_coordinates(poa)

# ==============================================================================
# GMSH model
# ==============================================================================

model = MeshModel.from_mesh(plate, targetlength=300)

model.heal()
model.refine_mesh()
model.generate_mesh(2)

# ==============================================================================
# COMPAS mesh
# ==============================================================================

compas_mesh = model.mesh_to_compas()
lengths = [compas_mesh.edge_length(edge) for edge in compas_mesh.edges()]

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ShellSection(material=mat, t=50)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec, implementation='shelldkgt')
mdl.add_part(prt)


# Set boundary conditions in the corners
for vertex in plate.vertices_where({'vertex_degree': 2}):
    location = plate.vertex_coordinates(vertex)
    mdl.add_pin_bc(nodes=prt.find_nodes_around_point(location, distance=10))


# Initialize a problem
prb = mdl.add_problem(name='SLS')
# Initialize a step
stp = prb.add_static_step()
# Create a load combination
stp.combination = LoadCombination.SLS()
# Add the loads
pt = prt.find_closest_nodes_to_point(poa_coordinates, distance=10)
stp.add_node_pattern(nodes=pt,
                      z=-1*units.kN,
                      load_case="LL")
stp.add_node_pattern(nodes=pt,
                      z=-(10*units.kN).to_base_units().magnitude,
                      load_case="DL")
stp.add_gravity_load_pattern([prt], g=9.81*units("m/s**2"), load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D', 'SF'])
stp.add_output(fout)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field

results_summary = {
    "Total volume ": mdl.volume,
    "Total weight ": mdl.volume*78/10**6*9.81/10,
    "Total vertical reaction": prb.get_total_reaction(stp)[0].z,
    "Reactions check": True if abs(round(prb.get_total_reaction(stp)[0].z-mdl.volume*23.5/10**6*9.81/10,0)) < 10 else False,
}
pprint(results_summary)

# Show Results
prb.show_displacements_contour(stp, scale_results=10, component=None, show_bcs=0.5)


