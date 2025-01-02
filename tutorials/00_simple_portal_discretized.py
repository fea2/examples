# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
from compas.geometry import Point

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import ElasticIsotropic, BeamElement, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.units import units
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])

# === Step 1: Define the Units System ===
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

compas_fea2.POINT_OVERLAP = False
# === Step 2: Create a Model and a Deformable Part ===
# Initialize the main finite element model
mdl = Model(name="discetized_portal")

# Create a deformable part that will contain nodes and elements
prt = DeformablePart(name="my_part")

# === Step 3: Define Material Properties ===
# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=30 * units("GPa"),       # Young's modulus (30 GPa)
    v=0.2,                     # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3")  # Density (2400 kg/m³)
)

# === Step 4: Create Nodes (Geometric Points) ===
# Define the nodes of the structure with their spatial coordinates
p1 = Point(x=0, y=0, z=0)       # Bottom-left corner of the structure
p2 = Point(x=0, y=0, z=3000)    # Top-left corner of the structure
p3 = Point(x=5000, y=0, z=3000) # Top-right corner of the structure
p4 = Point(x=5000, y=0, z=0)    # Bottom-right corner of the structure

# === Step 5: Define Cross-Sections ===
# Define rectangular cross-sections for columns and beams
sec_column = RectangularSection(
    w=20 * units.cm, h=30 * units.cm, material=mat  # Width = 20 cm, Height = 30 cm
)
sec_beam = RectangularSection(
    w=20 * units.cm, h=50 * units.cm, material=mat  # Width = 20 cm, Height = 50 cm
)

# === Step 6: Create Beam Elements ===
# Define column elements connecting vertical nodes
d = 10
nodes_column_1 = prt.add_nodes([Node([0, 0, i*(p2.z - p1.z)/d]) for i in range(d+1)])
nodes_column_2 = prt.add_nodes([Node([p4.x, 0, i*(p3.z - p4.z)/d]) for i in range(d+1)])
nodes_beam = prt.add_nodes([Node([i*(p3.x - p2.x)/d, 0, p2.z]) for i in range(d+1)])

for c, _ in enumerate(nodes_column_1[:-1]):
    prt.add_element(BeamElement(nodes=(nodes_column_1[c], nodes_column_1[c+1]),
                                                 section=sec_column,
                                                 frame=[0, 1, 0]))
for c, _ in enumerate(nodes_column_2[:-1]):
    prt.add_element(BeamElement(nodes=(nodes_column_2[c], nodes_column_2[c+1]),
                                                 section=sec_column,
                                                 frame=[0, 1, 0]))
for c, _ in enumerate(nodes_beam[:-1]):
    prt.add_element(BeamElement(nodes=(nodes_beam[c], nodes_beam[c+1]),
                                                 section=sec_beam,
                                                 frame=[0, 1, 0]))



# === Step 8: Add the Part to the Model ===
# Add the deformable part to the finite element model
mdl.add_part(part=prt)

# === Step 9: Summarize and Visualize the Model ===
# Print a summary of the model's components (nodes, elements, materials, etc.)
mdl.summary()

# Display a 3D visualization of the model

from compas.geometry import Plane
mdl.add_fix_bc(nodes=prt.find_nodes_on_plane(Plane.worldXY()))

# DEFINE THE PROBLEM
# define a step
stp = StaticStep()
stp.combination = LoadCombination.ULS()

stp.add_node_pattern(nodes=prt.find_closest_nodes_to_point([0, 0, 3000], distance=0.1),
                      x=1*units.kN,
                      load_case='LL')
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF'])
stp.add_output(fout)

# set-up the problem
prb = Problem('discretized_portal_Fx', mdl)
prb.add_step(stp)

mdl.add_problem(problem=prb)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

disp = prb.displacement_field 
react = prb.reaction_field
# print(disp.get_max_component(1, stp).magnitude)
print(react.get_max_result(1, stp).magnitude)
print(react.get_max_result(2, stp).magnitude)
print(react.get_max_result(3, stp).magnitude)

prb.show_reactions(stp, scale_results=2, show_bcs=0.5)
# prb.show_deformed(scale_results=1000, show_original=0.1, show_bcs=0.1)
# prb.show_displacements_contour(stp, scale_results=10, component=None, show_bcs=0.5)
