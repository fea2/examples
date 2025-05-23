# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
import compas_fea2
from compas_fea2.model import Model, Part, Node
from compas_fea2.model import ElasticIsotropic, BeamElement, RectangularSection
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
)
from compas_fea2.results import DisplacementFieldResults, SectionForcesFieldResults
from compas_fea2.units import units

compas_fea2.set_backend("compas_fea2_castem")

compas_fea2.POINT_OVERLAP = False

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# === Step 1: Define the Units System ===
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

# === Step 2: Create a Model and a Deformable Part ===
# Initialize the main finite element model
mdl = Model(name="portal")

# Create a deformable part that will contain nodes and elements
prt = Part(name="my_part")

# === Step 3: Define Material Properties ===
# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=30 * units("GPa"),  # Young's modulus (30 GPa)
    v=0.2,  # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3"),  # Density (2400 kg/m³)
)

# === Step 4: Create Nodes (Geometric Points) ===
# Define the nodes of the structure with their spatial coordinates
n1 = Node(xyz=[0, 0, 0])  # Bottom-left corner of the structure
n2 = Node(xyz=[0, 0, 3000])  # Top-left corner of the structure
n3 = Node(xyz=[5000, 0, 3000])  # Top-right corner of the structure
n4 = Node(xyz=[5000, 0, 0])  # Bottom-right corner of the structure

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
column_1 = BeamElement(
    nodes=[n1, n2],  # Nodes defining the column
    section=sec_column,  # Cross-section of the column
    frame=[0, 1, 0],  # Local frame direction (vertical)
)
column_2 = BeamElement(
    nodes=[n4, n3],  # Nodes defining the second column
    section=sec_column,
    frame=[0, 1, 0],
)

# Define a beam element connecting the top nodes
beam = BeamElement(
    nodes=[n2, n3],  # Nodes defining the beam
    section=sec_beam,  # Cross-section of the beam
    frame=[0, 1, 0],  # Local frame direction
)

# === Step 7: Add Elements to the Deformable Part ===
# Add the created elements (columns and beam) to the deformable part
prt.add_elements([column_1, column_2, beam])

# === Step 8: Add the Part to the Model ===
# Add the deformable part to the finite element model
mdl.add_part(part=prt)

# === Step 9: Summarize and Visualize the Model ===
# Print a summary of the model's components (nodes, elements, materials, etc.)
mdl.summary()

# Display a 3D visualization of the model
# mdl.show()

mdl.add_pin_bc(nodes=[n1, n4])

# DEFINE THE PROBLEM
prb = mdl.add_problem(problem=Problem(name="simple_portal_Fx"))
# define a step
stp = prb.add_step(StaticStep())
stp.combination = LoadCombination.ULS()
# Add a node pattern to apply a load on node n2
stp.add_uniform_node_load(nodes=[n2], x=1 * units.kN, load_case="LL")
stp.add_outputs([DisplacementFieldResults])

mdl.analyse_and_extract(problems=[prb], path=TEMP, output=True)
stp.show_deformed(scale_results=1000, show_original=0.1, show_bcs=0.1)
# beam.plot_stress_distribution(stp)
