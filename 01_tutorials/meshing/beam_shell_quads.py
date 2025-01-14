import os
from math import pi

from compas.datastructures import Mesh
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ShellSection, Steel
from compas_fea2.problem import (
    LoadCombination,
    DisplacementFieldOutput,
    ReactionFieldOutput,
    Stress2DFieldOutput,
)
from compas_fea2.units import units

units = units(system="SI_mm")

compas_fea2.set_backend("compas_fea2_opensees")

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1] + ["temp"])

lx = (3 * units.m).to_base_units().magnitude
ly = (100 * units.cm).to_base_units().magnitude
nx = 100
ny = 10
plate = Mesh.from_meshgrid(lx, nx, ly, ny).rotated(pi / 2, [1, 0, 0])

mdl = Model(name="plate_quads")
mat = Steel.S355()
sec = ShellSection(material=mat, t=30 * units.mm)
prt = mdl.add_part(
    DeformablePart.shell_from_compas_mesh(mesh=plate, section=sec, name="beam")
)

for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == 0:
        mdl.add_fix_bc(nodes=prt.find_nodes_around_point(location, distance=1))

prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_node_pattern(
    nodes=loaded_nodes,
    # y=-(2 / len(loaded_nodes)) * units.kN,
    z=-(20 / len(loaded_nodes)) * units.kN,
    load_case="LL",
)

# Ask for field outputs
# fout = FieldOutput(node_outputs=["U", "RF"], element_outputs=["S2D"])
stp.add_output(DisplacementFieldOutput())
stp.add_output(ReactionFieldOutput())
stp.add_output(Stress2DFieldOutput())

mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

cmap = ColorMap.from_color(Color.red(), rangetype="full")
# Show Results
# prt.sorted_nodes_by_displacement(stp)
# prb.show_principal_stress_vectors(stp, scale_results=0.01, show_bcs=0.05, show_loads=0.1)
# prb.show_deformed(scale_results=100000, show_nodes=True, fast=True, show_original=0.3, show_bcs=0.05, show_loads=0.1, opacity=0.8, original=0.25)
# prb.show_displacements_contour(
#     
# )
# # prb.show_stress_contour(stp, scale_results=0.5, show_bcs=0.05)

# prb.show_displacements(stp, fast=True, show_bcs=0.05, component=0, show_loads=0.1, show_contours=False, show_vectors=1000000, cmap=cmap )
prb.show_reactions(stp, show_bcs=0.05, show_vectors=0.1)
