import os
import compas
from compas.datastructures import Mesh
from compas.geometry import Scale, Plane

from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ShellSection, ElasticIsotropic
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")

# ==============================================================================
# Define the data files
# ==============================================================================
HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "00_data")
TEMP = os.path.join(HERE, "..", "temp")

cablemeshpath = os.path.join(DATA, "knitcandela", "CableMesh.json")
mesh: Mesh = compas.json_load(cablemeshpath)
# mesh.thickened((20 * units.mm).to_base_units().magnitude)

# Scale the geometry to mm
mesh.transform(Scale.from_factors([1000.0] * 3))

# ==============================================================================
# FEA2 Model
# ==============================================================================

# Initialize the Model object
mdl = Model(name="knitcandela")

# Use GMSH to discretize the geometry
print("Generating discretization...")
model = MeshModel.from_mesh(mesh, targetlength=500)
model.heal()
model.generate_mesh(2)
compas_mesh = model.mesh_to_compas()
print("Discretization complete!")

# Define mechanical properties and shell thickness
E = 10 * units.GPa
v = 0.2
rho = 1500 * units("kg/m**3")
t = 20 * units.mm

# Define material and section
mat = ElasticIsotropic(E=E, v=v, density=rho)
sec = ShellSection(t=t, material=mat)

# Define a deformable part using the mesh geometry and the mechanical properties
prt = Part.shell_from_compas_mesh(mesh=compas_mesh, section=sec)
mdl.add_part(prt)

# Fix the base
bottom_plane = Plane([0, 0, 0], [0, 0, 1])
fixed_nodes = prt.find_nodes_on_plane(bottom_plane, tol=10)
mdl.add_pin_bc(nodes=fixed_nodes)

# NOTE: if you run this, the results will not show
# mdl.show(fast=True, show_bcs=0.5)
# viewer = ModelViewer(mdl)
# viewer.show()

# ==============================================================================
# FEA2 Problem
# ==============================================================================
# Define the problem
prb = mdl.add_problem(name="VerticalLoad")

# Define the analysis step
stp = prb.add_static_step(system="SparseGeneral")
stp.combination = LoadCombination.SLS()

# Add a uniform load of 10 kN
loaded_nodes = prt.nodes
stp.add_node_pattern(
    nodes=loaded_nodes, load_case="LL", x=0, y=0, z=-10 * units.kN / len(loaded_nodes)
)

# Decide what information to save
stp.add_output(DisplacementFieldResults)

# Run the analysis and show results
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# # ==============================================================================
# # FEA2 Results
# # ==============================================================================

# NOTE: Currently the viewer can only show one plot each run
# NOTE: un-check the Model in the viewer to show the contour plot.
# stp.show_stress(stp, show_bcs=0.5, show_vectors=5000, show_contour=False, show_loads=10)
# stp.show_deformed(step=stp, scale_results=1000, show_original=0.5)
# stp.show_displacements(step=stp, fast=True, show_vectors=False, show_loads=10)
# stp.show_reactions(step=stp, show_bcs=0.5, show_vectors=3)

# viewer = ModelViewer(mdl)
# viewer.add_node_field_results(
#     stp.displacement_field, draw_cmap="viridis", draw_isolines=True
# )
# viewer.show()


viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
viewer.show()
