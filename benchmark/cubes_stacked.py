import os
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import FieldOutput, LoadCombination

from compas_fea2.units import units
units = units(system='SI_mm')

# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
step_file1 = os.path.join(DATA, "box.stp")
step_file2 = os.path.join(DATA, "box2.stp")


# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='cube')
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
box_bottom = DeformablePart.from_step_file(step_file=step_file1, **settings)
box_top = DeformablePart.from_step_file(step_file=step_file2, **settings)

mdl.add_part(box_bottom)
mdl.add_part(box_top)

# Set boundary conditions in the corners
restrained_nodes = list(filter(lambda n: n.z == 0 and (n.x == 0 or n.x==1000), box_bottom.nodes))
mdl.add_pin_bc(nodes=restrained_nodes)

prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = box_top.find_nodes_on_plane(box_top.top_plane) #list(filter(lambda n: n.z == 1000, prt.nodes))
stp.add_node_pattern(nodes=loaded_nodes,
                    z=-(1/len(loaded_nodes))*units.MN,
                    load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S3D'])
stp.add_output(fout)

mdl.add_problem(problem=prb)
mdl.summary()
prb.show(draw_bcs=0.1, draw_loads=1)

# Analyze and extracte results to SQLite database
mdl.analyse(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field


# Show Results
cmap = ColorMap.from_mpl('viridis')
prb.show_nodes_field_contour(disp, component=3, draw_reactions=0.01, draw_loads=0.01, draw_bcs=0.1, cmap=cmap)
# prb.show_nodes_field_vector(disp, component=3, scale_factor=1000, draw_bcs=0.01,  draw_loads=0.01)
# prb.show_deformed(scale_factor=1000, draw_bcs=0.01, draw_loads=0.01)
# prb.show_stress_contours(stress_type="von_mises_stress", draw_reactions=0.01, draw_loads=0.01, draw_bcs=0.1, cmap=cmap, bounds=[0, 0.5])
# prb.show_elements_field_vector(stress, vector_sf=100, draw_bcs=0.1)
