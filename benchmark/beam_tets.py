import os

from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput
from compas_fea2.results import NodeFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

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
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=(210*units.GPa).to_base_units().magnitude, 
                       v=0.2, 
                       density=(7800*units("kg/m**3")).to_base_units().magnitude)
sec = ShellSection(material=mat, t=100)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.shell_from_compas_mesh(mesh=plate, section=sec, name='beam', ndm=3, implementation="shelldkgq")
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == 0:
        mdl.add_fix_bc(nodes=prt.find_nodes_by_location(location, distance=1))

mdl.summary()
# mdl.show(draw_bcs=0.1)

# Initialize a step
stp = StaticStep()

# Add the load
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == lx:
        stp.add_node_load(nodes=prt.find_nodes_by_location(location, distance=1),
                          y=-(1*units.kN).to_base_units().magnitude)

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S', 'SF'])
stp.add_output(fout)

# Set-up the problem
prb = Problem('beam_shell')
prb.add_step(stp)
# prb.summary()
mdl.add_problem(problem=prb)

# Analyze and extracte results to SQLite database
# mdl.analyse(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

disp = NodeFieldResults(field_name='U', step=stp)
# print(disp.max.value)
# print(disp.min.value)


# # Show Results
# prb.show_nodes_field_contour('U', '2')
# prb.show_nodes_field_vector('U', vector_sf=500)
# prb.show_elements_field_vector('S', vector_sf=0.5)
# prb.show_deformed(scale_factor=500, draw_bcs=0.2, drawLoads=0.2, original=1.)
prb.show_deformed(scale_factor=100)

