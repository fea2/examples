import os

from compas.datastructures import Mesh
from compas.geometry import Scale, Plane

from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ShellSection, ElasticIsotropic
from compas_fea2.problem import LoadCombination, ModalAnalysis, LoadField, ConcentratedLoad
from compas_fea2.results import DisplacementFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# Set the backend implementation
compas_fea2.set_backend("compas_fea2_opensees")
# compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "00_data")
TEMP = os.path.join(HERE, "..", "temp")

mdl = Model(name="shell")

# Get the geometry from the obj file (and scale it to mm)
mesh = Mesh.from_obj(os.path.join(DATA, "shell", "tofea_f.obj"))
mesh.transform(Scale.from_factors([1000.0] * 3))

# Use GMSH to discretize the geometry
print("Generating discretization...")
model = MeshModel.from_mesh(mesh, targetlength=1000)
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
fixed_nodes = prt.find_nodes_on_plane(bottom_plane)
mdl.add_pin_bc(nodes=fixed_nodes)

# Define the problem
prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step(min_inc_size=0.1)
stp.combination = LoadCombination.SLS()
stp.add_uniform_node_load(nodes=prt.nodes, load_case="LL", x=0, y=0, z=-1 * units.kN)
# stp.add_gravity_load(parts=prt, g=9.81 * units("m/s**2"), z=-1., load_case="DL")

# Decide what information to save
stp.add_outputs([DisplacementFieldResults])

# stp_modal = ModalAnalysis(modes=3)
# prb.add_step(stp_modal)


# Run the analysis and show results
mdl.assign_keys()
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True, erase_data=True)
# # stp.show_displacements(step=stp, fast=True, show_vectors=0.5)
# stp_modal.show_mode_shape(
#     step=stp_modal,
#     mode=1,
#     scale_results=100,
#     show_original=0.3,
#     show_bcs=0.2,
#     show_vectors=100,
#     # show_contour=True,
# )

# Show Results
# Compas Viewer
stp.show_deformed(scale_results=1, show_original=0.3, show_bcs=0.5, show_loads=0.1)
# Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=1
)
viewer.show()
