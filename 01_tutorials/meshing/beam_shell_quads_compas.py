import os
from math import pi

from compas.datastructures import Mesh
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ShellSection, Steel
from compas_fea2.problem import (
    LoadCombination,
)
from compas_fea2.results import (
    DisplacementFieldResults,
    StressFieldResults,
    ReactionFieldResults,
)
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2] + ["temp"])

lx = (3 * units.m).to_base_units().magnitude
ly = (100 * units.cm).to_base_units().magnitude
nx = 100
ny = 30
plate = Mesh.from_meshgrid(lx, nx, ly, ny).rotated(pi / 2, [1, 0, 0])

mdl = Model(name="plate_quads")
mat = Steel.S355()
sec = ShellSection(material=mat, t=30 * units.mm)
prt = mdl.add_part(Part.shell_from_compas_mesh(mesh=plate, section=sec, name="beam"))

fixed_nodes = prt.nodes.subgroup(condition=lambda node: node.x == 0)
mdl.add_fix_bc(nodes=fixed_nodes)

prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

loaded_nodes = prt.nodes.subgroup(condition=lambda node: node.x == lx)
stp.add_uniform_node_load(
    nodes=loaded_nodes,
    # y=-(2 / len(loaded_nodes)) * units.kN,
    z=-(20 / len(loaded_nodes)) * units.kN,
    load_case="LL",
)
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])
# Ask for field outputs
# fout = FieldOutput(node_outputs=["U", "RF"], element_outputs=["S2D"])


mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True, erase_data=True)

# Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
# viewer.add_stress_tensors(stress, 1)
# viewer.add_principal_stress_vectors(stp.stress_field, 100)
viewer.show()


# Show Results
# prt.sorted_nodes_by_displacement(stp)
# prb.show_principal_stress_vectors(stp, scale_results=0.01, show_bcs=0.05, show_loads=0.1)
# prb.show_deformed(scale_results=100000, show_nodes=True, fast=True, show_original=0.3, show_bcs=0.05, show_loads=0.1, opacity=0.8, original=0.25)
# prb.show_displacements_contour(
#
# )
# prb.show_stress_contour(stp, scale_results=0.5, show_bcs=0.05)

# stp.show_displacements(fast=True, show_bcs=0.05, component=0, show_loads=0.1, show_contours=False, show_vectors=1000000, cmap=cmap )
# stp.show_reactions(show_bcs=0.05, show_vectors=0.1)
