import os
from math import pi
from compas.datastructures import Mesh

# from compas.utilities import geometric_key_xy
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import LoadCombination

from compas_fea2.results import StressFieldResults, DisplacementFieldResults

from compas_fea2.units import units

from compas_fea2_vedo.viewer import ModelViewer

units = units(system="SI_mm")

# compas_fea2.set_backend("compas_fea2_opensees")
# compas_fea2.set_backend("compas_fea2_calculix")
compas_fea2.set_backend("compas_fea2_castem")

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2] + ["temp"])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (1 * units.m).to_base_units().magnitude
ly = (15 * units.cm).to_base_units().magnitude
thk = (10 * units.mm).to_base_units().magnitude

plate = Mesh.from_polygons([[[0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0]]])
plate = plate.rotated(pi / 2, [1, 0, 0])
plate = plate.thickened(thk)
model = MeshModel.from_mesh(plate, targetlength=25)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name="beam_shell_tets")
# Define some properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))
sec = SolidSection(material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt = Part.from_gmsh(gmshModel=model, section=sec, name="beam")
mdl.add_part(prt)

# Set boundary conditions in the corners
fixed_nodes = list(filter(lambda n: n.x == 0, prt.nodes))
mdl.add_pin_bc(fixed_nodes)

mdl.summary()
# mdl.show(draw_bcs=0.1)

prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step(system="SparseGeneral", name="static")
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_uniform_node_load(
    nodes=loaded_nodes, z=-(1 / len(loaded_nodes)) * units.kN, load_case="LL"
)

# Ask for field outputs
stp.add_outputs([DisplacementFieldResults])
# prb.summary()

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, erase_data=True)
stress = stp.stress_field
disp = stp.displacement_field

# Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(disp, draw_cmap="viridis", draw_vectors=100)
# viewer.add_stress_tensors(stress, 1)
# viewer.add_principal_stress_vectors(stress, 100)
viewer.show()

# # print(stp.reaction_field.get_max_result("z").vector)
# fixed_nodes_top = [n for n in fixed_nodes if ly / 2 < n.z <= ly]
# fixed_nodes_bottom = [n for n in fixed_nodes if n.z <= ly / 2]
# # print(stp.reaction_field.compute_resultant(sub_set=fixed_nodes_bottom))
# # print(stp.reaction_field.compute_resultant(sub_set=fixed_nodes_top))
# resultant_F, resultant_M, loc = stp.reaction_field.compute_resultant(
#     sub_set=fixed_nodes
# )

# print(resultant_F)
# # Show Results
# stp.show_reactions(stp, show_vectors=0.1, show_bcs=0.05, show_contours=0.5)
# # stp.show_deformed(scale_results=1000, show_bcs=0.05, show_loads=0.1, show_original=0.25)
# # stp.show_displacements(
# #     show_vectors=1000, show_bcs=0.05, show_loads=0.1, show_contour=0.2
# # )
# # stp.show_stress(stp, scale_results=0.5, show_bcs=0.05, show_loads=0.1)
# # stp.show_stress(stp, show_bcs=0.05, show_vectors=100)
