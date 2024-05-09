import os

from compas.datastructures import Mesh

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
lx = (1*units.m).to_base_units().magnitude
ly = (30*units.cm).to_base_units().magnitude
nx = 20
ny = 4
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate_quads')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ShellSection(material=mat, t=30*units.mm)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.shell_from_compas_mesh(mesh=plate, section=sec, name='beam', ndm=3, implementation="ShellMIT4")
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == 0:
        mdl.add_fix_bc(nodes=prt.find_nodes_around_point(location, distance=1))

mdl.summary()

prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_node_pattern(nodes=loaded_nodes,
                  y=-(2/len(loaded_nodes))*units.kN,
                  load_case="LL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D'])
stp.add_output(fout)
# prb.summary()

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field
print(react.get_max_component(2, stp).magnitude)

# Show Results
prb.show_elements_field_vector(stress, vector_sf=1, draw_bcs=0.05, draw_loads=0.1)
prb.show_deformed(scale_factor=500, draw_bcs=0.05, draw_loads=0.1, opacity=0.8, original=0.25)
