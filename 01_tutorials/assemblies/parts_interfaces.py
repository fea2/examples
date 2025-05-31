import os
import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import SolidSection, Steel
from compas_fea2.units import units

from compas.geometry import Translation

# from compas_fea2.utilities.interfaces_numpy import mesh_mesh_interfaces

from compas_fea2_vedo.viewer import ModelViewer

units = units(system="SI_mm")

# ==============================================================================
# Define the data files
# ==============================================================================
HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "00_data", "solids")
TEMP = os.path.join(HERE, "..", "..", "temp")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")

mdl = Model(name="boxes")
mat_steel = Steel.S355(units=units)
sec = SolidSection(material=mat_steel)

prt = Part.from_step_file(
    step_file=os.path.join(DATA, "box.stp"),  # noqa: E501
    section=sec,
    meshsize_max=300,
)
mdl.add_part(prt)
print("Added part", 1)

w = prt.bounding_box.width
h = prt.bounding_box.height
for i in range(1, 4):
    mdl.copy_part(prt, Translation.from_vector([i * w, 0, 0]))
    print("Added part", i + 1)

print(set(mdl.graph.neighbors(prt)))
print(set(prt.graph.neighbors(list(prt.nodes)[10])))

part1 = list(mdl.parts)[0]
part2 = list(mdl.parts)[1]
part3 = list(mdl.parts)[2]


interfaces_dict = {}
for i, part in enumerate(mdl.parts):
    for j, neighbor in enumerate(mdl.parts):
        if i >= j:  # Avoid duplicate checks (A vs B is same as B vs A)
            continue

        mesh1 = part.bounding_box.to_mesh()
        mesh2 = neighbor.bounding_box.to_mesh()
        interfaces = mesh_mesh_interfaces(mesh1, mesh2)

        if interfaces:
            mdl.graph.add_edge(
                part, neighbor, relation="contact"
            )  # Add edge dynamically
            interfaces_dict[(part, neighbor)] = interfaces

print(interfaces_dict)
viewer = ModelViewer(mdl)
viewer.add_interfaces(interfaces=list(interfaces_dict.values()), color="teal")
viewer.show(show_parts=False)

