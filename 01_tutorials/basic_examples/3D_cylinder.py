# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os

from compas.geometry import Point, Line,  Cylinder

from compas_gmsh.models import ShapeModel

import compas_fea2
from compas_fea2.model import Model, Part, ElasticIsotropic, SolidSection
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2_vedo.viewer import ModelViewer
from compas_fea2.units import units


#--------------------------------------
# Backends, units and paths
#--------------------------------------
compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend("compas_fea2_sofistik")

units = units(system="SI")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")


# ==============================================================================
# Create a cylinder and meshing
# ==============================================================================
#geometry parameters
r = 1 #m
h = 4 #m

#compas.geometry objects
pb = Point(0,0,0)
ph = Point(0,0,h)
l_height = Line(pb, ph)

#definition of cylinder object in compas.geometry
#the centroid of the object is the origin of the frame O=(0, 0, 0,)
cylinder1 = Cylinder.from_line_and_radius(line=l_height, radius=r)

#mesh definition with compas_gmsh 
print('Beginning of gmsh meshing')
compasgmsh_model = ShapeModel(name='cylinder1')
compasgmsh_model.add_cylinder(cylinder1)
compasgmsh_model.options.mesh.lmax = h/20
compasgmsh_model.options.mesh.lmin = h/20
compasgmsh_model.generate_mesh(3)
print('End of gmsh meshing')

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================
# Initialize the model
mdl = Model(name="steel_shell")

# Define material properties
mat = ElasticIsotropic(E=300000 * units.Pa, 
                       v=0.3, 
                       density=2500 * units("kg/m**3"))

# Define the shell section
sec = SolidSection(material=mat)

# Create a deformable part from the mesh of gmsh
prt = Part.from_gmsh(gmshModel=compasgmsh_model, section=sec)
mdl.add_part(prt)

# Boundary conditions : the base is fixed
bottom_nodes = prt.nodes.subgroup(lambda n: n.z == 0)
mdl.add_fix_bc(nodes=bottom_nodes)

# Print model summary and visualize the model
# mdl.summary()
# mdl.show(show_bcs=0.001)

# ==============================================================================
# Define the problem
# ==============================================================================
# define step and combination
prb = mdl.add_problem(name="cylinder_traction")
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# define the loads
top_nodes = prt.nodes.subgroup(lambda n: n.z == h)
stp.add_uniform_node_load(nodes=top_nodes, z=1000 * units.N, load_case="LL")

# define the outputs
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdl.analyse_and_extract(
    problems=[prb], path=os.path.join(TEMP, prb.name), erase_data=True, verbose=True
)

# Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=100, show_original=0.1, show_bcs=0.5, show_loads=0.1)
#Vedo Viewer
viewer = ModelViewer(mdl)
viewer.add_node_field_results(stp.displacement_field, draw_cmap="viridis", draw_vectors=10)
# viewer.add_stress_tensors(stress, 1)
# viewer.add_principal_stress_vectors(stress, 100)
viewer.show()

# ==============================================================================
# Results verification
# ==============================================================================
# Force equilibrium

reaction_vector, applied_load = stp.check_force_equilibrium()
fz_tot = applied_load[2]

# Displacements
# z-axis (height)
disp = stp.displacement_field.get_result_at(list(top_nodes)[0])
max_z = stp.displacement_field.get_max_result("z")
print("Maximum displacement result in z (FEA) "+ str(max_z.z) +" mm")

area = 3.14 *r**2
eps_z = fz_tot/area/mat.E
z_th = eps_z*h
print("Theorical z-displacement value on the top "+ str(z_th) +" mm")


# r-axis (radius)
min_x = stp.displacement_field.get_min_result("x")
min_y = stp.displacement_field.get_min_result("y")
print("Minimum displacement result in x (FEA)  : "+ str(min_x.x)+" mm")

eps_r = -mat.v * eps_z
x_th = eps_r*r
print("Theorical x-displacement value : "+ str(x_th) +" mm")


