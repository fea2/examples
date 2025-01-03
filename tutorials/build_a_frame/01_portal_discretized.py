# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
import gmsh
from compas.geometry import Point, Line, Plane

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node, BeamElement
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.units import units
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, '..', '..', 'temp')

# === Step 1: Define the Units System ===
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

compas_fea2.POINT_OVERLAP = False

# === Step 2: Create a Model and a Deformable Part ===
# Initialize the main finite element model
mdl = Model(name="discretized_portal")

# Create a deformable part that will contain nodes and elements
prt = DeformablePart(name="my_part")

# === Step 3: Define Material Properties ===
# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=30 * units("GPa"),       # Young's modulus (30 GPa)
    v=0.2,                     # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3")  # Density (2400 kg/m³)
)

# === Step 4: Define Cross-Sections ===
# Define rectangular cross-sections for columns and beams
sec_column = RectangularSection(
    w=20 * units.cm, h=30 * units.cm, material=mat  # Width = 20 cm, Height = 30 cm
)
sec_beam = RectangularSection(
    w=20 * units.cm, h=50 * units.cm, material=mat  # Width = 20 cm, Height = 50 cm
)

# === Step 5: Create Geometry and Discretize with GMSH ===
# Initialize GMSH
gmsh.initialize()
gmsh.model.add("portal_frame")

# Define points for columns and beam
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(0, 0, 3000)
p3 = gmsh.model.geo.addPoint(5000, 0, 3000)
p4 = gmsh.model.geo.addPoint(5000, 0, 0)

# Define lines for columns and beam
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p3, p4)
l3 = gmsh.model.geo.addLine(p2, p3)

# Synchronize GMSH model
gmsh.model.geo.synchronize()

# Define mesh size and generate mesh
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 300)
gmsh.model.mesh.generate(1)

# Extract mesh data
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_tags, element_nodes = gmsh.model.mesh.getElementsByType(1)

# Create nodes and elements for the deformable part
nodes = []
for i in range(len(node_tags)):
    nodes.append(Node([node_coords[3 * i], node_coords[3 * i + 1], node_coords[3 * i + 2]]))

elements = []
for i in range(len(element_tags)):
    elements.append((int(element_nodes[2 * i] - 1), int(element_nodes[2 * i + 1] - 1)))

# Add nodes and elements to the part
prt.add_nodes(nodes)
for element in elements:
    if nodes[element[0]].z == nodes[element[1]].z:
        sec = sec_beam
    else:
        sec = sec_column
    prt.add_element(BeamElement(nodes=(nodes[element[0]], nodes[element[1]]), section=sec, frame=[0, 1, 0]))

# Finalize GMSH
gmsh.finalize()

# Add the part to the model
mdl.add_part(part=prt)

# === Step 6: Set Boundary Conditions ===
mdl.add_fix_bc(nodes=prt.find_nodes_on_plane(Plane.worldXY()))

# === Step 7: Define the Problem ===
# Define a static step and load combination
stp = StaticStep()
stp.combination = LoadCombination.ULS()

# Add a load at the top-left corner of the structure
stp.add_node_pattern(nodes=prt.find_closest_nodes_to_point([0, 0, 3000], distance=0.1), x=1 * units.kN, load_case='LL')

# Define field outputs
fout = FieldOutput(node_outputs=['U', 'RF'], element_outputs=['SF'])
stp.add_output(fout)

# Set up the problem
prb = Problem('discretized_portal_Fx', mdl)
prb.add_step(stp)

# Add the problem to the model
mdl.add_problem(problem=prb)

# === Step 8: Run the Analysis and Show Results ===
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Get displacement and reaction fields
disp = prb.displacement_field 
react = prb.reaction_field

# Print reaction results
print("Max reaction force in X direction [N]: ", react.get_max_result(1, stp).magnitude)
print("Max reaction force in Y direction [N]: ", react.get_max_result(2, stp).magnitude)
print("Max reaction force in Z direction [N]: ", react.get_max_result(3, stp).magnitude)

# Show reactions
prb.show_reactions(stp, scale_results=2, show_bcs=0.5)
