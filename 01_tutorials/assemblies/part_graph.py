import os
import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import SolidSection, Steel
from compas_fea2.units import units

from compas.geometry import Translation

units = units(system="SI_mm")

# ==============================================================================
# Define the data files
# ==============================================================================
HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "00_data", "solids")
TEMP = os.path.join(HERE, "..", "..", "temp")


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
