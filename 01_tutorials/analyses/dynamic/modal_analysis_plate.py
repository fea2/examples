"""
Tutorial 03: Modal Analysis of a Plate

This tutorial demonstrates how to perform a modal analysis on a simple plate
model using COMPAS FEA2.
We will define a plate made of a mesh, apply boundary conditions, and run a
modal analysis to extract and visualize the mode shapes.

Steps:
1. Define the plate geometry and material properties.
2. Create a deformable part from the mesh.
3. Set boundary conditions at both ends of the plate.
4. Define a modal analysis problem.
5. Run the analysis and visualize the mode shapes.
"""

import os

from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import Steel, SolidSection
from compas_fea2.problem import ModalAnalysis
from compas_fea2.units import units

from compas_fea2_vedo.viewer import ModelViewer

units = units(system="SI_mm")

# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "..", "temp")

# ==============================================================================
# Step 1: Define the plate geometry
# ==============================================================================
lx = (3 * units.m).to_base_units().magnitude
ly = (1 * units.m).to_base_units().magnitude
nx = 10
ny = 3
plate = Mesh.from_meshgrid(lx, nx, ly, ny)
plate = plate.thickened((1 * units.cm).to_base_units().magnitude)

# ==============================================================================
# Step 2: Create a deformable part from the mesh
# ==============================================================================
model = MeshModel.from_mesh(plate, targetlength=500)
solid_mesh = model.mesh_to_compas()

# Initialize the model
mdl = Model(name="modal_steel_plate")

# Define material properties
mat = Steel.S355()

# Define the solid section
sec = SolidSection(material=mat)

# Create a deformable part from the mesh
prt = Part.from_gmsh(gmshModel=model, section=sec)
mdl.add_part(prt)

# Assign mass to nodes for modal analysis
for n in prt.nodes:
    n.mass = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

# ==============================================================================
# Step 3: Set boundary conditions at both ends of the plate
# ==============================================================================
fixed_nodes = prt.nodes.subgroup(condition=lambda node: node.x == 0 or node.x == lx)
mdl.add_pin_bc(nodes=fixed_nodes)
# ==============================================================================
# Step 4: Define a modal analysis problem
# ==============================================================================
prb = mdl.add_problem(name="modal_analysis")
stp = prb.add_step(ModalAnalysis(modes=6))

# # ==============================================================================
# # Step 5: Run the analysis and visualize the mode shapes
# # ==============================================================================
prb.analyse_and_extract(path=TEMP, erase_data='armageddon')

# viewer = ModelViewer(mdl, shape=(2, 2))
# viewer.add_mode_shapes(list(stp.shapes)[:3], sf=1000)
# viewer.show()
