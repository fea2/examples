import os
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection, ZeroLengthSpringConnector, SpringSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')
# compas_fea2.VERBOSE = True
# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
box_1 = os.path.join(DATA, "box.stp")
box_2 = os.path.join(DATA, "box2.stp")


# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='stacked_boxes')
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


mdl = Model(name='stacked_boxes')
prt1 = DeformablePart.from_step_file(step_file=box_1, **settings)
prt2 = DeformablePart.from_step_file(step_file=box_2, **settings)

mdl.add_parts([prt1, prt2])

# Set boundary conditions in the corners
restrained_nodes = list(filter(lambda n: n.z == 0 and (n.x == 0 or n.x==1000), prt1.nodes))
mdl.add_pin_bc(nodes=restrained_nodes)

# Add connectors between boxes
top_1 = prt1.find_nodes_on_plane(prt1.top_plane)
bottom_2 =  prt2.find_nodes_on_plane(prt2.bottom_plane)
spring_section = SpringSection(axial=100, lateral=10, rotational=0)
for nt in top_1:
    nb = prt2.find_nodes_around_point(nt.point, distance=1, single=True)
    mdl.add_connector(ZeroLengthSpringConnector(nodes=[nb, nt], section=spring_section, directions=[0, 0, 1], failure={"c": 0}))

prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = prt1.find_nodes_on_plane(prt1.top_plane) #list(filter(lambda n: n.z == 1000, prt.nodes))
stp.add_node_pattern(nodes=loaded_nodes,
                    y=(1/len(loaded_nodes))*units.MN,
                    load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U'],
                #    element_outputs=['S3D']
                   )
stp.add_output(fout)

mdl.add_problem(problem=prb)
mdl.summary()
# mdl.show(draw_bcs=0.1, draw_loads=1)

# Analyze and extracte results to SQLite database
# mdl.analyse(problems=[prb], path=TEMP, verbose=True)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

stress = prb.stress_field

# prb.show_elements_field_vector(stress)
prb.show_deformed(scale_factor=1000)