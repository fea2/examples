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
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, SolidSection
from compas_fea2.problem import (
    LoadCombination,
)
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

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
mdl = Model(name="steel_plate")

# Define material properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))

# Define the solid section
sec = SolidSection(material=mat)

# Create a deformable part from the mesh
prt = Part.from_step_file(
    "/Users/francesco/code/fea2/examples/00_data/solids/box.stp",
    section=sec,
    meshsize_max=100,
)
mdl.add_part(prt)
mdl.show()
