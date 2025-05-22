"""
Tutorial 00: Beam Line Model

This tutorial demonstrates how to create a simple beam line model using COMPAS FEA2.
We will define a beam made of multiple line segments, apply boundary conditions,
add loads, and run a static analysis to extract and visualize the results.

Steps:
1. Define the beam geometry and material properties.
2. Create a deformable part from the line segments.
3. Set boundary conditions at the start of the beam.
4. Define a static problem and add loads.
5. Run the analysis and visualize the results.
"""

import os

from compas.geometry import Line

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_castem")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# ==============================================================================
# Create a beam line model
# ==============================================================================
# Define the length and number of segments
length = (0.1 * units.m).to_base_units().magnitude
n = 50

# Create lines representing the beam segments
lines = [Line([i * length, 0, 0], [(i + 1) * length, 0, 0]) for i in range(n - 1)]

# Initialize the model
mdl = Model(name="beam_line")

# Define material properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))

# Define the beam section
sec = ISection.HEA160(mat)

# Create a deformable part from the lines
prt = Part.from_compas_lines(lines, section=sec, name="beam")
mdl.add_part(prt)

# Set boundary conditions at the start of the beam
mdl.add_fix_bc(nodes=prt.nodes.subgroup(condition=lambda node: node.x == 0))

# Print model summary
mdl.summary()

# ==============================================================================
# Define the problem
# ==============================================================================
prb = mdl.add_problem(name="beam_line_Fz")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add a uniform distributed load
stp.add_uniform_node_load(
    nodes=prt.nodes.subgroup(condition=lambda node: node.x == (n - 1) * length),
    z=-1 * units.kN,
    load_case="LL",
)

# Define field outputs
stp.add_output(DisplacementFieldResults)
# mdl.show()
# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(
    problems=[prb], path=os.path.join(TEMP, prb.name), erase_data=True
)

# Get displacement and reaction fields
disp = stp.displacement_field


# Print displacement results
print("Displacement vector with min z component: ", disp.get_min_result("z").vector)
print("Resultant displacement: ", disp.compute_resultant())

# Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
viewer.show()
