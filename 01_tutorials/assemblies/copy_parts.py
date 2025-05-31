import os
import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import SolidSection, Steel
from compas_fea2.problem import Problem, LoadCombination
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2.units import units

from compas.geometry import Translation

from compas_fea2_vedo.viewer import ModelViewer

units = units(system="SI_mm")

# ==============================================================================
# Define the data files
# ==============================================================================
HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "00_data", "solids")
TEMP = os.path.join(HERE, "..", "..", "temp")

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

# ==============================================================================
# Define the model and add parts
# ==============================================================================
mdl = Model(name="boxes")
mat_steel = Steel.S355(units=units)
sec = SolidSection(material=mat_steel)

prt = Part.from_step_file(
    step_file=os.path.join(DATA, "box.stp"),  # noqa: E501
    section=sec,
    meshsize_max=100,
)
mdl.add_part(prt)
print("Added part", 1)

# copy the part multiple times to create an assembly
w = prt.bounding_box.width
h = prt.bounding_box.height
for i in range(1, 4):
    mdl.copy_part(prt, Translation.from_vector([i * w, 0, 0]))
    print("Added part", i + 1)

# Add a pin boundary condition to the bottom nodes
fixed_nodes = mdl.find_nodes_on_plane(mdl.bottom_plane, tol=10)
mdl.add_pin_bc(nodes=fixed_nodes)

# ==============================================================================
# Define the problem and analysis step
# ==============================================================================
prb = mdl.add_problem(Problem(name="distributed_load"))
stp_static = prb.add_static_step(name="static")
stp_static.combination = LoadCombination.SLS()

stp_static.add_gravity_load(parts=mdl.parts, g=9.81, name="gravity", load_case="DL")
for part in mdl.parts:
    loaded_nodes = part.find_nodes_on_plane(part.top_plane, tol=10)
    stp_static.add_uniform_node_load(
        name="uniform",
        nodes=loaded_nodes,
        z=-10 * units.kN / len(loaded_nodes),
        load_case="LL",
    )
stp_static.add_outputs([DisplacementFieldResults, ReactionFieldResults])

prb.analyse_and_extract(path=TEMP, erase_data=True, output=True)

# ==============================================================================
# Visualize the results
# ==============================================================================
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp_static.displacement_field,
    draw_vectors=10000,
    draw_cmap="viridis",
    draw_isolines=10,
)
viewer.show()
