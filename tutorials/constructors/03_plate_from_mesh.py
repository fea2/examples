"""
Tutorial 03: Plate from Mesh Model

This tutorial demonstrates how to create a simple plate model using COMPAS FEA2.
We will define a plate made of a mesh, apply boundary conditions,
add loads, and run a static analysis to extract and visualize the results.

Steps:
1. Define the plate geometry and material properties.
2. Create a deformable part from the mesh.
3. Set boundary conditions at both ends of the plate.
4. Define a static problem and add loads.
5. Run the analysis and visualize the results.
"""

import os

from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import LoadCombination, FieldOutput

from compas_fea2.units import units
units = units(system='SI_mm')

# Set the backend implementation
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, '..', '..', 'temp')

# ==============================================================================
# Create a plate model
# ==============================================================================
# Define the plate dimensions and mesh density
lx = (3 * units.m).to_base_units().magnitude
ly = (1 * units.m).to_base_units().magnitude
nx = 10
ny = 3
plate = Mesh.from_meshgrid(lx, nx, ly, ny)
plate = plate.thickened((1 * units.cm).to_base_units().magnitude)

# ==============================================================================
# GMSH model
# ==============================================================================
model = MeshModel.from_mesh(plate, targetlength=500)
solid_mesh = model.mesh_to_compas()

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================
# Initialize the model
mdl = Model(name='steel_plate')

# Define material properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))

# Define the solid section
sec = SolidSection(material=mat)

# Create a deformable part from the mesh
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec)
mdl.add_part(prt)

# Set boundary conditions at both ends of the plate
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == 0 or location[0] == lx:
        mdl.add_fix_bc(nodes=prt.find_nodes_around_point(location, distance=1))

# Print model summary
mdl.summary()

# ==============================================================================
# Define the problem
# ==============================================================================
prb = mdl.add_problem(name='mid_load')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add a load in the middle of the grid
loaded_nodes = []
for vertex in plate.vertices():
    location = plate.vertex_coordinates(vertex)
    if location[0] == lx / 2:
        loaded_nodes.extend(prt.find_nodes_around_point(location, distance=1))
stp.add_node_pattern(nodes=loaded_nodes, z=-10 * units.kN, load_case='LL')

# Define field outputs
fout = FieldOutput(node_outputs=['U', 'RF'], element_outputs=['S2D', 'SF'])
stp.add_output(fout)

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)

# Get displacement and reaction fields
disp = prb.displacement_field 
react = prb.reaction_field

# Print reaction and displacement results
print("Max reaction force [N]: ", react.get_max_result(3, stp).magnitude)
print("Min displacement [mm]: ", disp.get_min_component(3, stp))

# Show deformed shape
prb.show_deformed(scale_results=100, show_nodes=True, show_bcs=0.1)