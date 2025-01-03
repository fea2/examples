"""
Tutorial 00: Beam Line Model

This tutorial demonstrates how to create a simple beam line model using COMPAS FEA2.
We will define a beam made of multiple line segments, apply boundary conditions,
add loads, and run a modal analysis to extract and visualize the results.

Steps:
1. Define the beam geometry and material properties.
2. Create a deformable part from the line segments.
3. Set boundary conditions at the start of the beam.
4. Define a modal analysis problem.
5. Run the analysis and visualize the results.
"""

import os

from compas.geometry import Line

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import FieldOutput
from compas_fea2_opensees.problem.steps.perturbations import OpenseesModalAnalysis

from compas_fea2.units import units
units = units(system='SI_mm')

# Set the backend implementation
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, '..', '..', 'temp')

# ==============================================================================
# Create a beam line model
# ==============================================================================
# Define the length and number of segments
l = (0.1 * units.m).to_base_units().magnitude
n = 10

# Create lines representing the beam segments
lines = [Line([0, 0, i * l], [0, 0, (i + 1) * l]) for i in range(n)]

# Initialize the model
mdl = Model(name='beam_line')

# Define material properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))

# Define the beam section
sec = ISection(w=100, h=190, tw=2, tf=2, material=mat)

# Create a deformable part from the lines
prt = DeformablePart.from_compas_lines(lines, section=sec, name='beam', mass=10)
mdl.add_part(prt)

# Set boundary conditions at the start of the beam
mdl.add_fix_bc(nodes=prt.find_closest_nodes_to_point([0, 0, 0], distance=0.1))

# Apply mass at the top
top_node = prt.find_closest_nodes_to_point([0, 0, n * l], distance=0.1)[0]
top_node.mass = 100.0

# Print model summary
mdl.summary()

# ==============================================================================
# Define the problem
# ==============================================================================
prb = mdl.add_problem(name='beam_line_modal')
modal_step = OpenseesModalAnalysis(modes=5, name='modal_step')
prb.add_step(modal_step)

# Define field outputs
fout = FieldOutput(node_outputs=['U'], element_outputs=['SF'])
modal_step.add_output(fout)

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)

# # Get eigenvalues
# eigenvalues = prb.get_eigenvalues()
# print("Eigenvalues:", eigenvalues)

# # Show mode shapes
# prb.show_mode_shapes()