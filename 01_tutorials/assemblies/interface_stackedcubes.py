import os

from compas.geometry import Box, Vector, Translation

from compas_gmsh.models import ShapeModel

import compas_fea2
from compas_fea2.model import Model, Part, Interface, HardContactFrictionPenalty, SolidSection, ElasticIsotropic
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults, StressFieldResults

from compas_fea2_vedo import ModelViewer

from compas_fea2.units import units

#========================================
# Backend, paths & units
#========================================
compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_castem")

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "data")
TEMP = os.path.join(HERE, "..", "temp")

units = units(system="SI")

# ==============================================================================
# MODEL 
# ==============================================================================
w = 1 #m, width
h = 1 #m, height
d = 1 #m, depth

# Geometry and mesh with compas.geometry and compas_gmsh
box1 = Box.from_width_height_depth(w,h,d)
box1.translate(Vector(0,0,h/2))

model1 = ShapeModel(name='cube1')
model1.add_box(box1)
model1.options.mesh.lmin = 0.2
model1.generate_mesh(3)
compas_mesh1 = model1.mesh_to_compas()

# Define material and solid section (limestone)
mat = ElasticIsotropic(
E =  10* units("MPa"),
v = 0.2,
density = 2000 * units("kg/m**3"))

sec = SolidSection(mat)

# Creation of base part ()
prt = Part.from_boundary_mesh(
    compas_mesh1, section=sec, name="cube1"
)

#Compas_fea2 model
mdl = Model(name="stacked_cubes")
mdl.add_part(prt)

#cube dimensions
w = prt.bounding_box.width
h = prt.bounding_box.height

# Definition of the other part (cubes)
n_cubes = 3
for i in range(1,n_cubes):
    compas_mesh2=compas_mesh1.transformed(transformation=Translation.from_vector(vector=Vector(0,0,i*h)))
    mdl.add_part(Part.from_boundary_mesh(compas_mesh2, section=sec, name="cube"+str(i+1)))

# Behaviour law for the inteface
interaction = HardContactFrictionPenalty(mu=30, tolerance=1)

# Definition of interfaces
for i in range(1,len(mdl.parts)):
    part_below=mdl.find_part_by_name(name="cube"+str(i))
    part_above= mdl.find_part_by_name(name="cube"+str(i+1))

    slave_face = part_below.faces.subgroup(lambda n: n.centroid[2] == i*h)
    master_face = part_above.faces.subgroup(lambda n: n.centroid[2] == i*h)

    interface=Interface(master=master_face, slave=slave_face, behavior=interaction)
    mdl.add_interface(interface)

# Boundaries condition
fixed_nodes = prt.nodes.subgroup(lambda n: n.z == 0)
mdl.add_fix_bc(nodes=fixed_nodes)

# ==============================================================================
# PROBLEM 
# ==============================================================================
# static step and combination
prb = mdl.add_problem(name="SLS")
stp = prb.add_static_step(load_step=0.05, min_inc_size=0.1)
stp.combination = LoadCombination.SLS()

# loads : only the top cube is loaded with compression
stp.add_uniform_node_load(nodes=part_above.nodes.subgroup(lambda n: n.z == n_cubes*h), z=-10 * units.kN, load_case="LL")

# define the outputs
stp.add_outputs(
    [DisplacementFieldResults, StressFieldResults]
)

# ==============================================================================
# ANALYSIS
# ==============================================================================
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), output=True, max_increments=10, min_inc_size=0.01)

# ==============================================================================
# RESULTS AND VISUALIZATION
# ==============================================================================
# # Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=1000, show_original=0.3, show_loads=0.1)
#Compas Vedo
viewer = ModelViewer(model = mdl)
viewer.add_deformed_shape(step=stp,sf=1000)
viewer.show(show_parts=True, show_bcs=False)


