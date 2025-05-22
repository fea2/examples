import os

from random import choice, uniform
from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, ISection
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2_vedo.viewer import ModelViewer

from compas_fea2.units import units

units = units(system="SI_mm")

# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1] + ["temp"])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (10 * units.m).to_base_units().magnitude
ly = (10 * units.m).to_base_units().magnitude
nx = 8
ny = 5
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

inner_vertices = list(set(plate.vertices()) - set(plate.vertices_on_boundary()))

poa = choice(inner_vertices)
poa_coordinates = plate.vertex_coordinates(poa)

for v in inner_vertices:
    plate.vertex_attribute(
        v, "z", (uniform(-0.5, 0.5) * units.m).to_base_units().magnitude
    )

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name="grid")
# Define some properties
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))
sec = ISection(
    w=25 * units.cm, h=250 * units.mm, tw=20 * units.mm, tbf=4 * units.cm, ttf=4 * units.cm, material=mat
)

prt = Part.frame_from_compas_mesh(plate, section=sec)
mdl.add_part(prt)


# Set boundary conditions in the corners
for vertex in plate.vertices_where({"vertex_degree": 2}):
    location = plate.vertex_coordinates(vertex)
    mdl.add_pin_bc(nodes=prt.find_closest_nodes_to_point(location))

mdl.summary()

# Initialize a problem
prb = mdl.add_problem(name="SLS")
# Initialize a step
stp = prb.add_static_step()
# Create a load combination
stp.combination = LoadCombination.SLS()
# Add the loads
pt = prt.find_closest_nodes_to_point(poa_coordinates, 1)
stp.add_uniform_node_load(nodes=pt, z=-1 * units.kN, load_case="LL")
# stp.add_node_pattern(
#     nodes=pt, z=-(10 * units.kN).to_base_units().magnitude, load_case="DL"
# )
# stp.add_gravity_load_pattern([prt], g=9.81 * units("m/s**2"), load_case="DL")

# Ask for field outputs
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])

prb.summary()
# prb.show(show_bcs=1, draw_loads=0.1, opacity=1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = stp.displacement_field
react = stp.reaction_field
stress = stp.stress_field

results_summary = {
    "Total load": 1000.0,
    "Total vertical reaction": stp.get_total_reaction()[0].z,
    "Reactions check": (
        True
        if abs(
            round(
                stp.get_total_reaction()[0].z - 1000.0,
                0,
            )
        )
        == 0
        else False
    ),
}

print(results_summary)

# Show Results
#Compas Viewer
mdl.show()
stp.show_deformed(scale_results=0.5, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(
    stp.displacement_field, draw_cmap="viridis", draw_vectors=10000
)
viewer.show()
