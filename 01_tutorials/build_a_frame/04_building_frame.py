"""
Tutorial 02: 3D Frame Model

This tutorial demonstrates how to create a 3D frame model using COMPAS FEA2.
We will define a frame made of multiple beam elements, apply boundary conditions,
add loads, and run a static analysis to extract and visualize the results.

Steps:
1. Define the frame geometry and material properties.
2. Create a deformable part from the beam elements.
3. Set boundary conditions at the base of the frame.
4. Define a static problem and add loads.
5. Run the analysis and visualize the results.
"""

# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
import gmsh
from compas.geometry import Plane

import compas_fea2
from compas_fea2.model import Model, Part, Node, BeamElement
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
)
from compas_fea2.results import (
    DisplacementFieldResults,
    SectionForcesFieldResults,
    ReactionFieldResults,
)
from compas_fea2.units import units

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# === Step 1: Define the Units System ===
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

compas_fea2.POINT_OVERLAP = False

# === Step 2: Define Geometry Variables ===
# Define the dimensions of the frame
lx = 6000  # Length in X direction
ly = 4000  # Length in Y direction
lz = 3200  # Length in Z direction
nx = 2  # Number of divisions in X direction
ny = 3  # Number of divisions in Y direction
nz = 7  # Number of divisions in Z direction

# Define the target length for mesh discretization
target_length = 500

# === Step 3: Create a Model and a Deformable Part ===
# Initialize the main finite element model
mdl = Model(name="3d_frame")

# Create a deformable part that will contain nodes and elements
prt = Part(name="my_part")

# === Step 4: Define Material Properties ===
# Define an elastic isotropic material (e.g., concrete or steel)
mat = ElasticIsotropic(
    E=30 * units("GPa"),  # Young's modulus (30 GPa)
    v=0.2,  # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3"),  # Density (2400 kg/m³)
)

# === Step 5: Define Cross-Sections ===
# Define rectangular cross-sections for columns and beams
sec_column = RectangularSection(
    w=20 * units.cm, h=30 * units.cm, material=mat  # Width = 20 cm, Height = 30 cm
)
sec_beam = RectangularSection(
    w=20 * units.cm, h=50 * units.cm, material=mat  # Width = 20 cm, Height = 50 cm
)

# === Step 6: Create Geometry and Discretize with GMSH ===
# Initialize GMSH
gmsh.initialize()
gmsh.model.add("3d_frame")

# Define points for columns and beams
# Points are created in a 3D grid based on the number of divisions and dimensions
points = []
for i in range(nx + 1):
    for j in range(ny + 1):
        for k in range(nz + 1):
            points.append(gmsh.model.geo.addPoint(i * lx, j * ly, k * lz))

# Define lines for columns and beams
# Lines are created between the points to form the frame structure
lines = []
for i in range(nx + 1):
    for j in range(ny + 1):
        for k in range(nz):
            lines.append(
                gmsh.model.geo.addLine(
                    points[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k],
                    points[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k + 1],
                )
            )
for i in range(nx + 1):
    for j in range(ny):
        for k in range(nz + 1):
            lines.append(
                gmsh.model.geo.addLine(
                    points[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k],
                    points[i * (ny + 1) * (nz + 1) + (j + 1) * (nz + 1) + k],
                )
            )
for i in range(nx):
    for j in range(ny + 1):
        for k in range(nz + 1):
            lines.append(
                gmsh.model.geo.addLine(
                    points[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k],
                    points[(i + 1) * (ny + 1) * (nz + 1) + j * (nz + 1) + k],
                )
            )

# Synchronize GMSH model
gmsh.model.geo.synchronize()

# Define mesh size and generate mesh
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), target_length)
gmsh.model.mesh.generate(1)

# Extract mesh data
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_tags, element_nodes = gmsh.model.mesh.getElementsByType(1)

# Create nodes and elements for the deformable part
nodes = []
for i in range(len(node_tags)):
    nodes.append(
        Node([node_coords[3 * i], node_coords[3 * i + 1], node_coords[3 * i + 2]])
    )

elements = []
for i in range(len(element_tags)):
    elements.append((int(element_nodes[2 * i] - 1), int(element_nodes[2 * i + 1] - 1)))

# Add nodes and elements to the part
prt.add_nodes(nodes)
for element in elements:
    if nodes[element[0]].z == nodes[element[1]].z:
        if nodes[element[0]].z == 0:
            # prt.remove_node(nodes[element[0]])
            continue
        sec = sec_beam
        if nodes[element[0]].y == nodes[element[1]].y:
            frame = [0, 1, 0]
        else:
            frame = [1, 0, 0]
    else:
        sec = sec_column
        frame = [0, 1, 0]
    prt.add_element(
        BeamElement(
            nodes=(nodes[element[0]], nodes[element[1]]), section=sec, frame=frame
        )
    )

# Finalize GMSH
gmsh.finalize()

# Add the part to the model
mdl.add_part(part=prt)

# === Step 7: Set Boundary Conditions ===
mdl.add_fix_bc(nodes=prt.find_nodes_on_plane(Plane.worldXY()))

# === Step 8: Define the Problem ===
# Define a static step and load combination
stp = StaticStep()
stp.combination = LoadCombination.ULS()

# Add a load at the top-left corner of the structure
stp.add_uniform_node_load(
    nodes=prt.find_nodes_on_plane(Plane([0, 0, nz * lz], [0, 0, 1])),
    # x=1 * units.kN,
    # y=1 * units.kN,
    z=-1 * units.kN,
    load_case="LL",
)

# Define field outputs
fout = [DisplacementFieldResults, ReactionFieldResults]
stp.add_outputs(fout)

# Set up the problem
prb = Problem("3d_frame_Fx")
prb.add_step(stp)

# Add the problem to the model
mdl.add_problem(problem=prb)

# === Step 9: Run the Analysis and Show Results ===
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Get displacement and reaction fields
disp = stp.displacement_field
react = stp.reaction_field

# # Print reaction results
# print("Max reaction force in X direction [N]: ", react.get_min_result(1, stp).magnitude)
# print("Max reaction force in Y direction [N]: ", react.get_min_result(2, stp).magnitude)
# print("Max reaction force in Z direction [N]: ", react.get_max_result(3, stp).magnitude)

# Show reactions
# stp.show_displacements(stp, fast=True, show_bcs=0.5, show_loads=1, show_vectors=False)
stp.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=10)
stp.show_reactions(stp, show_vectors=0.05, show_bcs=0.05, show_contours=0.5)
