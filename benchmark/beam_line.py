import os

from compas.datastructures import Mesh
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (1*units.m).to_base_units().magnitude
ly = (20*units.cm).to_base_units().magnitude
nx = 20
ny = 4
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='beam_line')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = RectangularSection(w=10*units.cm, h=10*units.cm, material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.frame_from_compas_mesh(mesh=plate, section=sec, name='beam')
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == 0:
        mdl.add_fix_bc(nodes=prt.find_nodes_around_point(location, distance=1))

mdl.summary()
# mdl.show(draw_bcs=0.1)

prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes=[]
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == lx:
        loaded_nodes.append(*prt.find_nodes_around_point(location, distance=1))
stp.add_node_pattern(nodes=loaded_nodes, y=-1*units.kN, load_case='LL')


# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S', 'SF'])
stp.add_output(fout)

prb.show(draw_loads=0.1, draw_bcs=0.1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field
print(react.get_max_component(2, stp).magnitude)

prb.show_deformed(scale_factor=1000)