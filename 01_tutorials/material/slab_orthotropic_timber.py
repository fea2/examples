# IN THIS EXAMPLE, TWO WOODEN SHELL SLABS ARE MODELIZED
# SO AS TO SHOW THE IMPACT OF AN ORTHOTROPIC ELASTIC MATERIAL
# (TIMBER IS MODELIZED WITH SUCH A LAW), THESE WOODEN SHELL
# ARE SIMILAR EXCEPT FOR THE DIRECTION FOR THE GRAIN.
# THE SHELLS ARE IN CANTILEVER ALONG THE X-AXIS WITH A UNIFORM 
# LINE LOAD AT THE END OF THE SHELL.


# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os

from compas.datastructures import Mesh
from compas.geometry import Point, Line

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ShellSection, Timber
from compas_fea2.problem import LoadCombination
from compas_fea2.results import DisplacementFieldResults
from compas_fea2.units import units


#--------------------------------------
# Backend
#--------------------------------------
compas_fea2.set_backend("compas_fea2_abaqus")
# compas_fea2.set_backend("compas_fea2_sofistik")

units = units(system="SI_mm")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")


# ==============================================================================
# Create a shell model
# ==============================================================================
# Define the plate dimensions and mesh density
lx = (3 * units.m).to_base_units().magnitude
ly = (1 * units.m).to_base_units().magnitude
nx = 30
ny = 10
plate = Mesh.from_meshgrid(lx, nx, ly, ny)
thk = (10 * units.mm).to_base_units().magnitude

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================
# Initialize the model
mdlx = Model(name="slabx") #grain along x-axis
mdly = Model(name="slaby") #grain along y-axis

# Define material properties
mat = Timber.C24()

# Define the shells section
sec = ShellSection(t=thk, material=mat)

# Definition of the two parts, which longitudinal axis (grain axis) differs
prtx = Part.shell_from_compas_mesh(mesh=plate, section=sec, frame=[1,0,0])
mdlx.add_part(prtx)

prty = Part.shell_from_compas_mesh(mesh=plate, section=sec, frame=[0,1,0])
mdly.add_part(prty)


# The shells are fixed at x=0
fixed_nodes = prtx.nodes.subgroup(condition=lambda node: node.x == 0)
mdlx.add_fix_bc(nodes=fixed_nodes)

fixed_nodes = prty.nodes.subgroup(condition=lambda node: node.x == 0)
mdly.add_fix_bc(nodes=fixed_nodes)


# ==============================================================================
# Define the problem
# ==============================================================================

prbx = mdlx.add_problem(name="orthotropic")
stpx = prbx.add_static_step()
stpx.combination = LoadCombination.SLS()

prby = mdly.add_problem(name="mid_load")
stpy = prby.add_static_step()
stpy.combination = LoadCombination.SLS()

# Add a line load at the end of the slabs
loaded_nodes = prtx.nodes.subgroup(condition=lambda node: node.x == lx)
stpx.add_uniform_node_load(
    nodes=loaded_nodes, z=-1 * units.N / len(loaded_nodes), load_case="LL"
)

loaded_nodes = prty.nodes.subgroup(condition=lambda node: node.x == lx)
stpy.add_uniform_node_load(
    nodes=loaded_nodes, z=-1 * units.N / len(loaded_nodes), load_case="LL"
)

# Define field outputs
stpx.add_output(DisplacementFieldResults)
stpy.add_output(DisplacementFieldResults)

# ==============================================================================
# Run the analysis and show results
# ==============================================================================
# Analyze and extract results to SQLite database
mdlx.analyse_and_extract(
    problems=[prbx], path=str(os.path.join(TEMP, prbx.name)), output=True
)

mdly.analyse_and_extract(
    problems=[prby], path=str(os.path.join(TEMP, prby.name)), output=True
)

# Results
min_z_x = stpx.displacement_field.get_min_result("z").z
min_z_y = stpy.displacement_field.get_min_result("z").z

print(f"""In this example, the orthotropic material behaviour might be observed through deflection.
When the grain (which is the stiffer direction) is along the direction of the cantilever,
deflection is lower than when the grain is perpendicular to this direction.

Here the cantilever direction is along x-axis.
      
Grain along x-axis deflection: {abs(min_z_x)}mm
Grain along y-axis deflection: {abs(min_z_y)}mm""")
# # Show deformed shape
# #Compas Viewer
stpy.show_deformed(scale_results=1000, show_bcs=0.5, show_loads=0.1)

# #Plot deflection mid-width
line = Line(Point(0,ly/2,0), Point(lx,ly/2,0))
stpx.plot_deflection_along_line(line, n_divide=100)