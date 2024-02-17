import os

from compas.datastructures import Mesh
from compas.utilities import geometric_key_xy
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput
from compas_fea2.results import NodeFieldResults, ElementFieldResults

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

plate = Mesh.from_polygons([[[0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0]]])
model = MeshModel.from_mesh(plate, targetlength=50)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ShellSection(material=mat, t=30*units.mm)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec, name='beam', ndm=3, implementation='shelldkgt')
prt._discretized_boundary_mesh = model.mesh_to_compas()
prt._boundary_mesh = plate
prt.bounding_box
mdl.add_part(prt)

# Set boundary conditions in the corners
for node in prt.nodes:
    if node.x == 0:
        mdl.add_fix_bc(nodes=[node])

mdl.summary()
# mdl.show(draw_bcs=0.1)

# Initialize a step
stp = StaticStep()

# Add the load
loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_node_load(nodes=loaded_nodes,
                  y=-(2/len(loaded_nodes))*units.kN)

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D'])
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

stress_field = ElementFieldResults('S2D', stp)

for r in stress_field.results:
    print(r.von_mises_stress)
# print(stress_field.results[10].smin)

# # Show Results
# prb.show_nodes_field_contour('U', '2')
# prb.show_nodes_field_vector('U', vector_sf=500)
prb.show_elements_field_vector('S2D', vector_sf=1, draw_bcs=0.05, draw_loads=0.1)
prb.show_deformed(scale_factor=500, draw_bcs=0.05, draw_loads=0.1, opacity=0.8, original=0.25)
