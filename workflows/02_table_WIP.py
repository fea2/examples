import os

from random import choice
from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel
from compas.colors import ColorMap, Color
from compas.geometry import Point, Line

import compas_fea2
from compas_fea2.model import Model, DeformablePart, BeamElement, Node
from compas_fea2.model import ElasticIsotropic, ShellSection, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (1.5*units.m).to_base_units().magnitude
ly = (3*units.m).to_base_units().magnitude
nx = 5
ny = 10
plate = Mesh.from_meshgrid(lx, nx, ly, ny)

# ==============================================================================
# Select random internal vertex for load application
# ==============================================================================

poa = choice(list(set(plate.vertices()) - set(plate.vertices_on_boundary())))
poa_coordinates = plate.vertex_coordinates(poa)

# ==============================================================================
# GMSH model
# ==============================================================================

model = MeshModel.from_mesh(plate, targetlength=300)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# model.optimize_mesh(niter=100)
# model.recombine_mesh()

# ==============================================================================
# COMPAS mesh
# ==============================================================================

compas_mesh = model.mesh_to_compas()
lengths = [compas_mesh.edge_length(edge) for edge in compas_mesh.edges()]

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ShellSection(material=mat, t=50)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec, implementation='shelldkgt')
mdl.add_part(prt)

sec_beam = RectangularSection(w=50*units.mm, h=100*units.mm, material=mat)

# Convert the gmsh model in a compas_fea2 Part
h = (700*units.mm).to_base_units().magnitude
nh = 7

leg_A = [Line(Point(0, 0, -h+c*h/nh), Point(0, 0, -h+(c+1)*h/nh)) for c in range(nh)]
leg_B = [Line(Point(lx, 0, -h+c*h/nh), Point(lx, 0, -h+(c+1)*h/nh)) for c in range(nh)]
leg_C = [Line(Point(lx, ly, -h+c*h/nh), Point(lx, ly, -h+(c+1)*h/nh)) for c in range(nh)]
leg_D = [Line(Point(0, ly, -h+c*h/nh), Point(0, ly, -h+(c+1)*h/nh)) for c in range(nh)]
legs = [leg_A, leg_B, leg_C, leg_D]
for lines in legs:
    for line in lines:
        #FIXME change tolerance
        nodes = [prt.find_nodes_around_point(list(p), 1, single=True) or Node(list(p)) for p in list(line)]
        prt.add_nodes(nodes)
        element = BeamElement(nodes=nodes, section=sec_beam)
        prt.add_element(element)

A = mdl.find_nodes_around_point(Point(0, 0, -h), distance=1)
B = mdl.find_nodes_around_point(Point(lx, 0, -h), distance=1)
C = mdl.find_nodes_around_point(Point(lx, ly, -h), distance=1)
D = mdl.find_nodes_around_point(Point(0, ly, -h), distance=1)

# Set boundary conditions in the corners
mdl.add_pin_bc(nodes=A+B+C+D)

mdl.summary()
mdl.show(draw_bcs=0.3)

# Initialize a problem
prb = mdl.add_problem(name='SLS')
# Initialize a step
stp = prb.add_static_step()
# Create a load combination
stp.combination = LoadCombination.SLS()
# Add the loads
pt = prt.find_closest_nodes_to_point(poa_coordinates, distance=10)
stp.add_node_pattern(nodes=pt,
                      z=-1*units.kN,
                      load_case="LL")
# stp.add_node_pattern(nodes=pt,
#                       z=-(10*units.kN).to_base_units().magnitude,
#                       load_case="DL")
stp.add_gravity_load_pattern([prt], g=9.81*units("m/s**2"), load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D', 'SF'])
stp.add_output(fout)

# #TODO check how to apply loads in steps in OpenSees
# # stp2 = prb.add_static_step()
# # stp2.add_node_pattern(nodes=pt,
# #                       z=-(10*units.kN).to_base_units().magnitude,
# #                       load_case="DL")
# # # Create a load combination
# # stp2.combination = LoadCombination.SLS()

# # prb.summary()
# prb.show(draw_loads=0.1)

# Analyze and extracte results to SQLite database
# mdl.analyse(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field

print("\nDisplacement vectors with max and min component along 2:")
for r in disp.get_limits_component(2, step=stp):
    print(r.vector)

print("\nAbsolute max and min displacement vectors :")
for r in disp.get_limits_absolute(step=stp):
    print(r.vector)

# # # cmap = ColorMap.from_color(Color.red(), rangetype='light')
# # cmap = ColorMap.from_mpl('viridis')

# # # # Show Results
prb.show_deformed(scale_factor=100, draw_bcs=0.1, draw_loads=0.1)
# # prb.show_nodes_field_contour(disp, component=3, draw_reactions=0.2, draw_loads=0.5, draw_bcs=0.1, cmap=cmap)
# # # prb.show_stress_contours(stress_type="von_mises_stress", side="top", draw_reactions=0.02, draw_loads=0.05, draw_bcs=0.5, cmap=cmap)
# # # prb.show_elements_field_vector(stress, vector_sf=10, draw_bcs=1)
# # # prb.show_nodes_field_vector(field_name='U', scale_factor=500, draw_bcs=0.5,  draw_loads=0.1)
# # # prb.show_deformed(scale_factor=10, draw_bcs=0.5, draw_loads=0.1)


