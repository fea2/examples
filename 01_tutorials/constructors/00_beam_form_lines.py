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
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import LoadCombination, DisplacementFieldOutput

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")

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
sec = ISection(w=100, h=190, tw=2, tf=2, material=mat)

# Create a deformable part from the lines
prt = DeformablePart.from_compas_lines(lines, section=sec, name="beam")
mdl.add_part(prt)

# Set boundary conditions at the start of the beam
mdl.add_fix_bc(nodes=prt.find_closest_nodes_to_point([0, 0, 0], 1))

# Print model summary
mdl.summary()

# ==============================================================================
# Define the problem
# ==============================================================================
prb = mdl.add_problem(name="beam_line_Fz")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add a load at the end of the beam
stp.add_node_pattern(
    nodes=prt.find_closest_nodes_to_point([(n - 1) * length, 0, 0], 1),
    z=-1 * units.kN,
    load_case="LL",
)

# Define field outputs
fout = DisplacementFieldOutput()
stp.add_output(fout)

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)

# # Get displacement and reaction fields
# disp = prb.displacement_field

# # Print displacement results
# print("Displacement vector with min z component: ", disp.get_min_result(3, stp).vector)
# print(
#     "Min z component of the displacement field [mm]: ", disp.get_min_component(3, stp)
# )

# Show deformed shape
prb.show_displacements(
    stp, show_bcs=0.1, component=2, show_vectors=10, show_contours=False, show_loads=0.2
)
# prb.show_deformed(scale_results=10, show_original=0.1, show_bcs=0.1)
