# Import necessary classes from compas_fea2 for creating the model, materials, and elements
import os
import compas_fea2

from compas.geometry import Point, Line

from compas_fea2.model import Model, Part, ElasticIsotropic, RectangularSection
from compas_fea2.problem import Problem, StaticStep, LoadCombination
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2.units import units

#--------------------------------------
# Backend
#--------------------------------------
compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend("compas_fea2_sofistik")

HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "..", "temp")

#--------------------------------------
# Units
#--------------------------------------
# Define the unit system to be used (SI with millimeters)
units = units(system="SI_mm")  # SI units with length in millimeters

#--------------------------------------
# Model
#--------------------------------------
# Initialize the main finite element model
mdl = Model(name="simplebeam")

# Geometry with compas.geometry objects
# The geometry being defined through compas.geometry, the units pint methods can't be used
# The users must care to the consistency of the units
lx = 5000 #mm 

p1 = Point(x = 0, y=0, z=0)
p2 = Point(x = lx , y=0, z=0) 
l1 = Line(p1, p2)

# Material and section

mat = ElasticIsotropic(
    E=30 * units("GPa"),  # Young's modulus (30 GPa)
    v=0.2,  # Poisson's ratio (dimensionless)
    density=2400 * units("kg/m**3"),  # Density (2400 kg/m³)
)

sec = RectangularSection(
    w= 12 * units.cm, #width
    h=10 * units.cm, # height
    material= mat 
)

# Meshing and creation of part with specific method
prt = Part.from_compas_lines_discretized(lines=[l1], targetlength=10, element_model='BeamElement', section=sec, frame=[0,1,0])
mdl.add_part(prt)

# Meshing and creation of part with specific method
mdl.add_pin_bc(nodes=prt.find_closest_nodes_to_point(p1))
mdl.add_pin_bc(nodes=prt.find_closest_nodes_to_point(p2))

# Visualize the geometry 
# mdl.show()

#--------------------------------------
# Problem
#--------------------------------------

# define the problem
prb = mdl.add_problem(problem=Problem(name="simplebeam"))
# define a step
stp = prb.add_step(StaticStep())
stp.combination = LoadCombination.SLS()
# Add a uniform load
stp.add_uniform_node_load(nodes= mdl.nodes, z=-1 * units.kN, load_case="LL")
# define the outputs
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])

#--------------------------------------
# Analysis
#--------------------------------------

mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

#--------------------------------------
# Visualization
#--------------------------------------
disp = stp.displacement_field

stp.show_deformed(scale_results=10, show_original=0.1, show_bcs=0.1)

# stp.show_displacements()

# Print reaction results
react = stp.reaction_field
print("Max/Min reaction forces in Z direction [N]: ", 
      react.get_limits_component('z')[0].magnitude, "/",
      react.get_limits_component('z')[1].magnitude)
#--------------------------------------
# VERIFICATION OF RESULTS
#--------------------------------------

reaction_vector, applied_load = stp.check_force_equilibrium()

#result from FEA

max_z = stp.displacement_field.get_min_result("z")
print("Deflection result from FEA model : " +str(max_z.z) + " mm")

#theory
load_linear = applied_load[2]/lx
fl = 5/384 * load_linear * lx**4 / (mat.E * sec.Ixx)
print("Theorical deflection result : " +str(fl) + " mm")

