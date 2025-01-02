import os

from random import choice
from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel
from compas.colors import ColorMap, Color
from compas.geometry import Point, Line

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults


from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (10*units.m).to_base_units().magnitude
ly = (10*units.m).to_base_units().magnitude
nx = 8
ny = 5
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

poa = choice(list(set(plate.vertices()) - set(plate.vertices_on_boundary())))
poa_coordinates = plate.vertex_coordinates(poa)



# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='grid')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ISection(w=25*units.cm, h=250*units.mm, tw=20*units.mm, tf=4*units.cm, material=mat)

prt = DeformablePart.frame_from_compas_mesh(plate, section=sec)
mdl.add_part(prt)


# Set boundary conditions in the corners
for vertex in plate.vertices_where({'vertex_degree': 2}):
    location = plate.vertex_coordinates(vertex)
    mdl.add_pin_bc(nodes=prt.find_nodes_around_point(location, distance=10))

mdl.summary()

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

prb.summary()
# prb.show(show_bcs=1, draw_loads=0.1, opacity=1)

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


# Show Results
prb.show_displacements_contour(stp, scale_results=0.5, component=None, show_bcs=0.5)

