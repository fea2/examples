"""
Tutorial 01: Beam Grid from Mesh

This tutorial demonstrates how to create a simple grid model using COMPAS FEA2.
We will define a grid made of multiple beam segments, apply boundary conditions,
add loads, and run a static analysis to extract and visualize the results.

Steps:
1. Define the grid geometry and material properties.
2. Create a deformable part from the mesh.
3. Set boundary conditions at both ends of the beam.
4. Define a static problem and add loads.
5. Run the analysis and visualize the results.
"""

import os

from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# ==============================================================================
# Create a beam lines model
# ==============================================================================
# Define the plate dimensions and mesh density
lx = (3 * units.m).to_base_units().magnitude
ly = (1 * units.m).to_base_units().magnitude
nx = 10
ny = 3
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================
# Initialize the model
mdl = Model(name="steel_grid")

# Define material properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))

# Define the beams section
sec = ISection.HEA180(mat)

# Create a deformable part from the mesh
prt = Part.frame_from_compas_mesh(mesh=plate, section=sec, name="grid")
mdl.add_part(prt)

# Set boundary conditions at both ends of the beam
fixed_nodes = prt.nodes.subgroup(condition=lambda node: node.x == 0 or node.x == lx)
mdl.add_fix_bc(nodes=fixed_nodes)

# Print model summary
mdl.summary()

# ==============================================================================
# Define the problem
# ==============================================================================
prb = mdl.add_problem(name="mid_load")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add a load in the middle of the grid
loaded_nodes = prt.nodes.subgroup(condition=lambda node: node.x == lx / 2)
stp.add_node_pattern(nodes=loaded_nodes, z=-1 * units.kN, load_case="LL")

# Define field outputs

stp.add_output(DisplacementFieldResults)

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)

# # Get displacement and reaction fields
# disp = prb.displacement_fields
# react = prb.reaction_field

# # Print reaction and displacement results
# print("Max reaction force [N]: ", react.get_max_result(2, stp).magnitude)
# print("Min displacement [mm]: ", disp.get_min_result(2, stp).magnitude)

# Show deformed shape
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
viewer.show()
