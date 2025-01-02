import os
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

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


mdl = Model(name='box')
prt = DeformablePart.from_step_file(step_file=step_file, **settings)

print(prt._boundary_mesh)
mdl.add_part(prt)

# Set boundary conditions in the corners
restrained_nodes = list(filter(lambda n: n.z == 0 and (n.x == 0 or n.x==1000), prt.nodes))
mdl.add_pin_bc(nodes=restrained_nodes)

prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = prt.find_nodes_on_plane(prt.top_plane) #list(filter(lambda n: n.z == 1000, prt.nodes))
stp.add_node_pattern(nodes=loaded_nodes,
                    z=-(1/len(loaded_nodes))*units.MN,
                    load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S3D'])
stp.add_output(fout)

mdl.add_problem(problem=prb)
mdl.summary()
# mdl.show(draw_bcs=0.1, draw_loads=1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field


# Show Results
prb.show_principal_stress_vectors(stp, scale_results=20, show_bcs=0.05)
# prb.show_deformed(scale_results=100, show_bcs=0.05, show_loads=0.1, opacity=0.8, original=0.25)
# prb.show_displacements_contour(stp, scale_results=0.5, show_bcs=0.05, component=0)
# prb.show_stress_contour(stp, scale_results=0.5, show_bcs=0.05)
# prb.show_reactions(stp, scale_results=0.005, show_bcs=0.05)