import os
import compas_fea2
from random import choice
from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel
from compas.geometry import Point, Line
from compas_fea2.model import (Model, DeformablePart, BeamElement, Node,
                                ElasticIsotropic, ShellSection, RectangularSection)
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.units import units

# Initialize units and backend
units = units(system='SI_mm')
compas_fea2.set_backend('compas_fea2_opensees')

# Paths
HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, 'temp')

# Constants
LX = 1.5 * units.m
LY = 3 * units.m
NX = 5
NY = 10
H = 700 * units.mm
NH = 7

# ==============================================================================
# Plate Mesh
# ==============================================================================
plate = Mesh.from_meshgrid(LX.to_base_units().magnitude, NX, 
                           LY.to_base_units().magnitude, NY)

# Random vertex for load application
poa = choice(list(set(plate.vertices()) - set(plate.vertices_on_boundary())))
poa_coordinates = plate.vertex_coordinates(poa)

# ==============================================================================
# GMSH Model
# ==============================================================================
model = MeshModel.from_mesh(plate, targetlength=300)
model.heal()
model.refine_mesh()
model.generate_mesh(2)

compas_mesh = model.mesh_to_compas()

# ==============================================================================
# COMPAS_FEA2 Model
# ==============================================================================
# Initialize model
mdl = Model(name='plate')

# Material and sections
mat = ElasticIsotropic(E=210 * units.GPa, v=0.2, density=7800 * units("kg/m**3"))
shell_sec = ShellSection(material=mat, t=50)
beam_sec = RectangularSection(w=50 * units.mm, h=100 * units.mm, material=mat)

# Add deformable part from GMSH model
prt = DeformablePart.from_gmsh(gmshModel=model, section=shell_sec, implementation='shelldkgt')
mdl.add_part(prt)

# ==============================================================================
# Generate Supporting Legs
# ==============================================================================
def generate_legs(base_coords, height, sections):
    """Generate vertical legs for the structure at specified base points."""
    legs = []
    for base in base_coords:
        for i in range(sections):
            start = Point(base[0], base[1], -height + i * height / sections)
            end = Point(base[0], base[1], -height + (i + 1) * height / sections)
            legs.append(Line(start, end))
    return legs

# Base points for the legs
base_points = [
    (0, 0),  # Bottom-left corner
    (LX.to_base_units().magnitude, 0),  # Bottom-right corner
    (LX.to_base_units().magnitude, LY.to_base_units().magnitude),  # Top-right corner
    (0, LY.to_base_units().magnitude),  # Top-left corner
]
legs = generate_legs(base_points, H.to_base_units().magnitude, NH)

# Add all legs to the model
for leg in legs:
    nodes = []
    for p in [leg.start, leg.end]:
        # Find or create nodes for the leg
        node = prt.find_nodes_around_point([p.x, p.y, p.z], 1, single=True)
        if not node:
            node = Node([p.x, p.y, p.z])
            prt.add_node(node)
        nodes.append(node)
    # Add the beam element for each leg
    prt.add_element(BeamElement(nodes=nodes, section=beam_sec, frame=[0, 1, 0]))

# Boundary Conditions
corner_points = [Point(*pt, -H.to_base_units().magnitude) for pt in base_points]
bc_nodes = [mdl.find_nodes_around_point(pt, distance=1) for pt in corner_points]
mdl.add_pin_bc(nodes=[node for group in bc_nodes for node in group])

# ==============================================================================
# Problem Definition
# ==============================================================================
prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()

# Load Combination
stp.combination = LoadCombination.SLS()

# Add loads
load_nodes = prt.find_closest_nodes_to_point(poa_coordinates, distance=10)
stp.add_node_pattern(nodes=load_nodes, z=-1 * units.kN, load_case="LL")
stp.add_gravity_load_pattern([prt], g=9.81 * units("m/s**2"), load_case="DL")

# Field Outputs
fout = FieldOutput(node_outputs=['U', 'RF'], element_outputs=['S2D', 'SF'])
stp.add_output(fout)

# ==============================================================================
# Analysis
# ==============================================================================
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

# Results
disp = prb.displacement_field
print("\nMax and Min Displacement Component (Z):")
for r in disp.get_limits_component(2, step=stp):
    print(r.vector)

print("\nAbsolute Max and Min Displacements:")
for r in disp.get_limits_absolute(step=stp):
    print(r.vector)

# Visualization
prb.show_displacements_contour(stp, scale_results=100, component=None, show_bcs=0.5)