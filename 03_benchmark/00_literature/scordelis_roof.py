import os
import gmsh
from math import radians, cos, sin
from compas.geometry import Plane
import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node, ShellElement
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import (
    Problem,
    StaticStep,
    LoadCombination,
    DisplacementFieldOutput,
)
from compas_fea2.units import units

# --------------------------------------------------------------------------
# Initialize COMPAS FEA2
# --------------------------------------------------------------------------
compas_fea2.set_backend("compas_fea2_opensees")
units = units(system="SI_mm")
HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2] + ["temp"])

# --------------------------------------------------------------------------
# GMSH Geometry
# --------------------------------------------------------------------------
# Initialize GMSH
gmsh.initialize()
gmsh.model.add("Scordelis-Lo Roof")

# Parameters
R = 25000  # Radius in mm
L = 50000  # Length in mm
theta = 40  # Angle in degrees
t = 250  # Thickness in mm
n_angle = 20  # Mesh divisions along the curve
n_span = 20  # Mesh divisions along the span
q = -900 * 47.8803  # Uniform load in N/m^2 (converted to N/mm^2)

# Create points on the curved surface
angle_divisions = [radians(theta) * i / (n_angle - 1) for i in range(n_angle)]
span_divisions = [L * i / (n_span - 1) for i in range(n_span)]

nodes = []
for span in span_divisions:
    for angle in angle_divisions:
        x = R * sin(angle)
        y = span
        z = R * (1 - cos(angle))
        nodes.append(gmsh.model.geo.addPoint(x, y, z))

# Create lines and define surfaces
element_surfaces = []
for i in range(n_span - 1):
    for j in range(n_angle - 1):
        n1 = nodes[i * n_angle + j]
        n2 = nodes[i * n_angle + j + 1]
        n3 = nodes[(i + 1) * n_angle + j + 1]
        n4 = nodes[(i + 1) * n_angle + j]

        # Create lines for the element
        l1 = gmsh.model.geo.addLine(n1, n2)
        l2 = gmsh.model.geo.addLine(n2, n3)
        l3 = gmsh.model.geo.addLine(n3, n4)
        l4 = gmsh.model.geo.addLine(n4, n1)

        # Add curve loop and surface
        curve_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        surface = gmsh.model.geo.addPlaneSurface([curve_loop])
        element_surfaces.append(surface)

# Synchronize the geometry
gmsh.model.geo.synchronize()

# Mesh generation
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 2000)
gmsh.model.mesh.generate(2)

# --------------------------------------------------------------------------
# Convert GMSH Mesh to COMPAS FEA2
# --------------------------------------------------------------------------
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_tags, element_node_tags = gmsh.model.mesh.getElementsByType(
    2
)  # Triangular elements


def nodes_from_gmsh(node_tags, node_coords):
    nodes = {}
    for i, tag in enumerate(node_tags):
        x, y, z = node_coords[i * 3 : (i + 1) * 3]
        nodes[tag] = Node((x, y, z))
    return nodes


def elements_from_gmsh(element_tags, element_node_tags, nodes, section):
    elements = []
    for i in range(len(element_tags)):
        element_nodes = [nodes[n] for n in element_node_tags[i * 3 : (i + 1) * 3]]
        elements.append(ShellElement(nodes=element_nodes, section=section))
    return elements


# Convert nodes and elements
gmsh_nodes = nodes_from_gmsh(node_tags, node_coords)
material = ElasticIsotropic(name="Concrete", E=30e3, v=0.0, density=2500)
shell_section = ShellSection(name="ShellSection", t=t, material=material)

# Create model and part
mdl = Model(name="Scordelis-Lo Roof")
part = DeformablePart(name="RoofPart")
mdl.add_part(part)

part.add_nodes(list(gmsh_nodes.values()))
shell_elements = elements_from_gmsh(
    element_tags, element_node_tags, gmsh_nodes, shell_section
)
part.add_elements(shell_elements)

gmsh.finalize()

# --------------------------------------------------------------------------
# Boundary Conditions
# --------------------------------------------------------------------------
# Fix the curved edges
curved_edges = part.find_nodes_on_plane(
    Plane((0, 0, 0), (1, 0, 0))
) + part.find_nodes_on_plane(Plane((0, 0, 0), (-1, 0, 0)))
mdl.add_fix_bc(nodes=curved_edges)

# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------
problem = Problem(name="ScordelisLoProblem")
mdl.add_problem(problem)
step = problem.add_step(StaticStep(name="LoadStep"))
step.add_outputs([DisplacementFieldOutput()])

# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
# Apply uniform load to the surface
load_nodes = part.nodes
distributed_load = {node: (0, 0, q) for node in load_nodes}
step.combination = LoadCombination.ULS()

# Run analysis
mdl.analyse_and_extract(problems=[problem], path=TEMP, verbose=True)

# Show Results
# prb.show_reactions(stp, show_vectors=0.1, show_bcs=0.05, show_contours=0.5)
problem.show_deformed(
    scale_results=10000, show_bcs=1, show_original=0.25, show_loads=10
)
# prb.show_displacements(show_vectors=1000, show_bcs=0.05, show_loads=0.1, show_contour=0.2)
# prb.show_principal_stress_vectors(stp, scale_results=0.5, show_bcs=0.05, show_loads=0.1)
# prb.show_stress_contour(stp, scale_results=0.5, show_bcs=0.05)
