from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import os

from compas.geometry import Plane

import compas_fea2
from compas_fea2.model import Model, Part, Node
from compas_fea2.model import CircularSection, ElasticIsotropic, BeamElement
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
)
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")
compas_fea2.set_backend("compas_fea2_opensees")

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1] + ["temp"])

mdl = Model(name="ILOVEFEA2")
prt = Part(name="my_part")

mat = ElasticIsotropic(E=210 * units("GPa"), v=0.2, density=7800 * units("kg/m**3"))

# Define the word and font
word = "ILOVEFEA2"
font = FontProperties(
    family="sans-serif", size=(3000 * units.mm).to_base_units().magnitude
)

sec = CircularSection(r=10 * units.cm, material=mat)

# Get the path for the word
nodes_dict = {}
for c, letter in enumerate(word):
    text_path = TextPath(
        (c * (2000 * units.mm).to_base_units().magnitude, 0), letter, prop=font
    )
    vertices = text_path.vertices

    for start, end in zip(vertices, vertices[1:-1]):
        start_key = tuple(start)
        end_key = tuple(end)

        if start_key not in nodes_dict:
            nodes_dict[start_key] = prt.add_node(Node(xyz=[start[0], 0, start[1]]))
        if end_key not in nodes_dict:
            nodes_dict[end_key] = prt.add_node(Node(xyz=[end[0], 0, end[1]]))

        prt.add_element(
            BeamElement(
                nodes=[nodes_dict[start_key], nodes_dict[end_key]],
                section=sec,
                frame=[0, 1, 0],
            )
        )

mdl.add_part(part=prt)

fixed_nodes = prt.find_nodes_on_plane(plane=Plane.worldXY())

mdl.add_fix_bc(nodes=fixed_nodes)

# DEFINE THE PROBLEM
# define a step
stp = StaticStep()
stp.combination = LoadCombination.ULS()
stp.add_node_pattern(nodes=prt.nodes, y=1 * units.kN, load_case="LL")
stp.add_output(DisplacementFieldResults)

# set-up the problem
prb = Problem(name="fea2model")
prb.add_step(stp)

mdl.add_problem(problem=prb)

mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
viewer.show()
