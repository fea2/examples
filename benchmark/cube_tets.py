import os

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput
from compas_fea2.results import NodeFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
step_file = os.path.join(DATA, "box.stp")


# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='box')
# Define some properties
mat = ElasticIsotropic(E=30*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = SolidSection(material=mat)

# Convert the gmsh model in a compas_fea2 Part
settings = {       
    'target_mesh_size': 200,
    'mesh_size_at_vertices': None,
    'target_point_mesh_size': None,
    'meshsize_max': 300,
    'meshsize_min': 200,
    'rigid': False,
    "material": mat,
    "section": sec,
}


mdl = Model(name='box')
prt = DeformablePart.from_step_file(step_file=step_file, **settings)

mdl.add_part(prt)

# Set boundary conditions in the corners
restrained_nodes = list(filter(lambda n: n.z == 0 and (n.x == 0 or n.x==1000), prt.nodes))
mdl.add_pin_bc(nodes=restrained_nodes)

# Initialize a step
stp = StaticStep()

# Add the load
loaded_nodes = prt.find_nodes_on_plane(prt.top_plane) #list(filter(lambda n: n.z == 1000, prt.nodes))
stp.add_node_load(nodes=loaded_nodes,
                  z=-(1/len(loaded_nodes))*units.MN)

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S3D'])
stp.add_output(fout)

# Set-up the problem
prb = Problem('top_load')
prb.add_step(stp)
print(stp)
mdl.add_problem(problem=prb)
mdl.summary()
# mdl.show(draw_bcs=0.1, draw_loads=1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# disp = NodeFieldResults(field_name='U', step=stp)
# # print(disp.max.value)
# # print(disp.min.value)


# # Show Results
# prb.show_nodes_field_contour('U', '2')
# prb.show_nodes_field_vector('RF', vector_sf=0.001, draw_loads=0.002, draw_bcs=0.05)
prb.show_elements_field_vector('S3D', vector_sf=100, draw_bcs=0.05, draw_loads=0.002)
# prb.show_deformed(scale_factor=1000, draw_bcs=0.2, draw_loads=0.002, original=0.1)
# prb.show_deformed()

