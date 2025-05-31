import os
from math import pi
from compas.datastructures import Mesh
from random import choice

# from compas.utilities import geometric_key_xy
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import LoadCombination

from compas_fea2.results import (
    StressFieldResults,
    DisplacementFieldResults,
    ReactionFieldResults,
)

from compas_fea2.units import units

from compas_fea2_vedo.viewer import ModelViewer

units = units(system="SI_mm")

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2] + ["temp"])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (1 * units.m).to_base_units().magnitude
ly = (30 * units.cm).to_base_units().magnitude

plate = Mesh.from_polygons([[[0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0]]])
plate = plate.rotated(pi / 2, [1, 0, 0])
model = MeshModel.from_mesh(plate, targetlength=50)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name="beam_shell_tri")
# Define some properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))
sec = ShellSection(material=mat, t=30 * units.mm)

# Convert the gmsh model in a compas_fea2 Part
prt = Part.from_gmsh(gmshModel=model, section=sec, name="beam")
prt._discretized_boundary_mesh = model.mesh_to_compas()
prt._boundary_mesh = plate
prt.bounding_box
mdl.add_part(prt)

# Set boundary conditions in the corners
mdl.add_fix_bc(nodes=prt.nodes.subgroup(lambda n: n.x == 0))

mdl.summary()
mdl.show(draw_bcs=0.1)
# viewer = ModelViewer(mdl)
# viewer.show()

prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step(system="SparseGeneral")
stp.combination = LoadCombination.SLS()

# Add the load
loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_uniform_node_load(
    nodes=loaded_nodes, z=-(2 / len(loaded_nodes)) * units.kN, load_case="LL"
)

# Ask for field outputs
stp.add_outputs([DisplacementFieldResults])
# prb.summary()

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True, erase_data=True)
# print(react.get_max_result(2, stp).magnitude)

# Show Results
#Compas Viewer
stp.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_vectors=100000, draw_cmap="viridis"
)
# viewer.add_principal_stress_vectors(stp.stress_field, 1)
viewer.show()
