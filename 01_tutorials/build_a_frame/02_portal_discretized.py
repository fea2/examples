# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
import gmsh
from compas.geometry import Plane, Point, Line

import compas_fea2
from compas_fea2.model import Model, Part, BeamElement
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
)
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2.units import units

compas_fea2.set_backend("compas_fea2_castem")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# === Step 1: Define the Units System ===
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

compas_fea2.POINT_OVERLAP = False

# === Step 2: Create a Model and a Deformable Part ===
# Initialize the main finite element model
mdl = Model(name="discretized_portal")

# === Step 3: Define Material Properties ===
# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=30 * units("GPa"),  # Young's modulus (30 GPa)
    v=0.2,  # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3"),  # Density (2400 kg/m³)
)

# === Step 4: Define Cross-Sections ===
# Define rectangular cross-sections for columns and beams
sec_column = RectangularSection(
    w=20 * units.cm, h=30 * units.cm, material=mat  # Width = 20 cm, Height = 30 cm
)
sec_beam = RectangularSection(
    w=20 * units.cm, h=50 * units.cm, material=mat  # Width = 20 cm, Height = 50 cm
)

# === Step 5: Create Geometry and Discretize ====
#Create the geometry with compas.geometry

p1 = Point(0, 0, 0)
p2 = Point(0, 0, 3000)
p3 = Point(5000, 0, 3000)
p4 = Point(5000, 0, 0)

l1 = Line(p1, p2)
l2 = Line(p3, p4)
l3 = Line(p2, p3)

#
prt = Part.from_compas_lines_discretized(lines=[l1, l3, l2], targetlength=1000, element_class=BeamElement, section=sec_beam, frame=[0,1,0], name='portal')

# Add the part to the model
mdl.add_part(part=prt)

# === Step 6: Set Boundary Conditions ===
mdl.add_fix_bc(nodes=prt.find_nodes_on_plane(Plane.worldXY()))

# mdl.show()
# === Step 7: Define the Problem ===
# Set up the problem
prb = mdl.add_problem(problem=Problem(name="discretized_portal_Fx"))
stp = prb.add_step(StaticStep(name='horizontal_load'))
# Define a static step and load combination
stp.combination = LoadCombination.ULS()

# Add a load at the top-left corner of the structure
stp.add_uniform_node_load(
    nodes=prt.find_closest_nodes_to_point(point=p2, number_of_nodes=1),
    x=10 * units.kN,
    load_case="LL",
)

# Define field outputs
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])

# === Step 8: Run the Analysis and Show Results ===
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Print reaction results
react = stp.reaction_field
print("Max/Min reaction forces in X direction [N]: ", 
      react.get_limits_component('x')[0].magnitude, "/",
      react.get_limits_component('x')[1].magnitude)
print("Max/Min reaction forces in Y direction [N]: ", 
      react.get_limits_component('y')[0].magnitude, "/",
      react.get_limits_component('y')[1].magnitude)
print("Max/Min reaction forces in Y direction [N]: ", 
      react.get_limits_component('z')[0].magnitude, "/",
      react.get_limits_component('z')[1].magnitude)

# Show reactions
stp.show_deformed(scale_results=1000, show_original=0.1, show_bcs=0.1)
stp.show_reactions(stp, scale_results=2, show_bcs=0.5)