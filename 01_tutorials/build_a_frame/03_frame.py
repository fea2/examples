import os
import gmsh
from compas.geometry import Plane, cross_vectors, normalize_vector
import compas_fea2
from compas_fea2.model import Model, Part, Node, BeamElement, ShellElement
from compas_fea2.model import ElasticIsotropic, ShellSection, ISection, Steel
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
)

from compas_fea2.results import DisplacementFieldResults, StressFieldResults
from compas_fea2.units import units

# --------------------------------------------------------------------------
# 1. Initialize COMPAS FEA2 and Units
# --------------------------------------------------------------------------
# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_calculix")
# compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend('compas_fea2_sofistik')

units = units(system="SI_mm")
compas_fea2.POINT_OVERLAP = False

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

# --------------------------------------------------------------------------
# 2. Create Model and Deformable Part
# --------------------------------------------------------------------------
mdl = Model(name="table_with_shell")
prt = Part(name="table_part")
mdl.add_part(prt)
# --------------------------------------------------------------------------
# 3. Define Material and Sections
# --------------------------------------------------------------------------

steel = Steel.S355()
conc = ElasticIsotropic(
    E=30 * units("GPa"),  # Steel: 210 GPa
    v=0.17,  # Poisson's ratio
    density=2500 * units("kg/m**3"),  # Steel density
)

# Use ISection profiles
top_frame_sec = ISection.IPE200(steel)
leg_sec = ISection.HEA140(steel)
shell_section = ShellSection(
    t=100 * units.mm, material=conc
)  # Tabletop shell thickness

# --------------------------------------------------------------------------
# 4. Define Geometry in GMSH
# --------------------------------------------------------------------------
gmsh.initialize()
gmsh.model.add("table_3d")

# Dimensions
Lx, Ly, H = 4000, 5000, 2750

# Define coordinates for the top points
top_coords = [
    (0, 0, H),  # Point 1
    (Lx, 0, H),  # Point 2
    (Lx, Ly, H),  # Point 3
    (0, Ly, H),  # Point 4
]

# Create top points and store their tags
top_points = [gmsh.model.geo.addPoint(*coord) for coord in top_coords]

# Create bottom points using the same X, Y coordinates but Z=0
bottom_coords = [(x, y, 0) for x, y, _ in top_coords]
bottom_points = [gmsh.model.geo.addPoint(*coord) for coord in bottom_coords]

# Top frame lines
top_lines = [
    gmsh.model.geo.addLine(top_points[i], top_points[(i + 1) % 4]) for i in range(4)
]

# Leg lines
leg_lines = [gmsh.model.geo.addLine(top_points[i], bottom_points[i]) for i in range(4)]

# Create tabletop shell
top_loop = gmsh.model.geo.addCurveLoop(top_lines)
top_surface = gmsh.model.geo.addPlaneSurface([top_loop])

gmsh.model.geo.synchronize()

# Mesh settings
mesh_size = 100  # Approximate mesh size in mm
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
gmsh.model.mesh.generate(2)

# --------------------------------------------------------------------------
# 5. Convert GMSH Mesh to COMPAS FEA2 Elements
# --------------------------------------------------------------------------
# Extract nodes and elements
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
beam_element_tags, beam_element_nodes = gmsh.model.mesh.getElementsByType(1)  # Lines
shell_element_tags, shell_element_nodes = gmsh.model.mesh.getElementsByType(2)  # Shells

# Add nodes to the model
nodes = []
for i, tag in enumerate(node_tags):
    slice_idx = 3 * i
    slice_end = slice_idx + 3
    x, y, z = node_coords[slice_idx:slice_end]
    nodes.append(Node([x, y, z]))
prt.add_nodes(nodes)

# Add beam elements (legs and frame) with correct frame orientation
for i in range(len(beam_element_tags)):
    n1_id = int(beam_element_nodes[2 * i] - 1)
    n2_id = int(beam_element_nodes[2 * i + 1] - 1)
    n1, n2 = nodes[n1_id], nodes[n2_id]

    # Determine section and frame orientation
    axis_vector = normalize_vector([n2.x - n1.x, n2.y - n1.y, n2.z - n1.z])
    if abs(axis_vector[2]) > max(abs(axis_vector[0]), abs(axis_vector[1])):  # Vertical
        section = leg_sec
        frame = [1, 0, 0]  # Local Y-axis
    else:  # Horizontal
        section = top_frame_sec
        frame = normalize_vector(
            cross_vectors(axis_vector, [0, 0, 1])
        )  # Perpendicular in plane

    prt.add_element(BeamElement(nodes=(n1, n2), section=section, frame=frame))

# Add shell elements (tabletop)
for i in range(len(shell_element_tags)):
    elem_node_ids = [
        int(shell_element_nodes[3 * i] - 1),
        int(shell_element_nodes[3 * i + 1] - 1),
        int(shell_element_nodes[3 * i + 2] - 1),
    ]
    elem_nodes = [nodes[n_id] for n_id in elem_node_ids]
    prt.add_element(ShellElement(nodes=elem_nodes, section=shell_section))

gmsh.finalize()

# --------------------------------------------------------------------------
# 6. Assign Boundary Conditions
# --------------------------------------------------------------------------
mdl.add_fix_bc(nodes=prt.find_nodes_on_plane(Plane.worldXY()))  # Fix bottom nodes

# mdl.show()

# --------------------------------------------------------------------------
# 7. Define Problem, Load, and Outputs
# --------------------------------------------------------------------------
prb = mdl.add_problem(Problem(name="table_problem"))
stp = prb.add_step(StaticStep(min_inc_size=0.01))
stp.combination = LoadCombination.ULS()

# Add a distributed load on the tabletop shell
load_nodes = prt.find_nodes_on_plane(Plane([0, 0, H], [0, 0, 1]))
stp.add_uniform_node_load(
    nodes=load_nodes, z=-3 * units.kN , load_case="LL"
)

# Add output fields
stp.add_outputs([DisplacementFieldResults])
# --------------------------------------------------------------------------
# 8. Run Analysis and Show Results
# --------------------------------------------------------------------------
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
# mdl.show()
# stp.show_stress(stp, show_bcs=0.05, show_vectors=1e2, plane="bottom")

stp.show_deformed(scale_results=1000, show_original=0.2, show_bcs=0.3, show_loads=10)
# prb.show_displacements(stp, show_bcs=0.3, show_loads=10, show_contour=True, show_vectors=100)
