import compas_fea2
import os

from compas.geometry import Frame, Point, Line

from compas_fea2.model import Model, RectangularSection, ElasticIsotropic, RigidLinkConnector, Part
from compas_fea2.problem import Problem, StaticStep, LoadCombination
from compas_fea2.results import DisplacementFieldResults, ReactionFieldResults
from compas_fea2.units import units


#========================================
# Backend, paths & units
#========================================
# Initialization of plugins and paths
compas_fea2.set_backend("compas_fea2_castem")
# compas_fea2.set_backend("compas_fea2_abaqus")
compas_fea2.POINT_OVERLAP = False
HERE = os.path.dirname(__file__)
TEMP = os.path.join(HERE, "..", "temp")

# Unit system
units = units(system="SI_mm")  # SI units with length in millimeters

# Define an elastic isotropic material (timber, C24)
mat = ElasticIsotropic(
    E=11 * units("GPa"),  # Young's modulus (11 GPa)
    v=0.35,  # Poisson's ratio (dimensionless)
    density=420 * units("kg/m**3"),  # Density (420 kg/m³)
)

#========================================
# MODEL
#========================================
# Define the model
mdl = Model(name="roof_timber")

# Geometry
p00 = Point(0, 0, 0)
p03 = Point(3000, 0, 0)
p06 = Point(6000, 0, 0)
p13 = Point(3000, 0, 1000)
p22 = Point(2000, 0, 2000)
p24 = Point(4000, 0, 2000)
p33 = Point(3000, 0, 3000)

# Lines
lines = {
    "Part0": Line(p00, p06), #Part 0 - tie beam
    "Part1": Line(p03, p33), #Part 1 -  king post
    "Part2": Line(p00, p33), #Part 2 - left principal rafter
    "Part3": Line(p33, p06), #Part 3 - right principal rafter
    "Part4": Line(p22, p13), #Part 4 - left strut
    "Part5": Line(p13, p24)  #Part 5 - right strut
}

# Section
sec = RectangularSection(
    w=10 * units.cm,
    h=30 * units.cm,
    material=mat
)

# Constructor for determining the frame of a beam
def frame_oriented(line, gamma=0):
    """
    Determine the local frame of a line according to its direction.
    The rotation according to the neutral-axis can be changed with the gamma parameters (radian).

    Parameters
    ----------
    line : :class: compas.geometry.Line
            Compas line defining the neutral axis of the beam.
    
    gamma : float (optional)
            Rotation of the frame along the longitudinal axis.
            If not indicated, gamma is considered null.

    """
    from math import atan, pi, asin
    from numpy import sign

    v = line.vector
    if v.y== 0 :
        alpha = -pi/2 *sign(v.x)
        if v.x==0:
            beta = pi/2
        else :
            beta = asin(v.z/v.length)
    else :
        alpha = -atan(v.x/v.y) 

        beta = asin(v.z/v.length)*v.y/abs(v.y)
    
    f = Frame.from_euler_angles(
        euler_angles=[alpha , beta, gamma],
        static=False,
        axes="zxy",
        point=line.start,
    )
    return f

# part definition from discretizing method
parts = []
for partname, line in lines.items() :
    prt = mdl.add_part(Part.from_compas_lines_discretized(lines=[line], targetlength=500, element_model='BeamElement', section=sec, frame=frame_oriented(line, gamma=0), name=partname))
    parts.append(prt)

# Definition of the connection between the parts
# list with data for connector definition between parts
connection_data= [
    #1st part   #2nd part    #common point  #description
    [0,         1,          p03],           #tie beam/king post
    [2,         0,          p00],           #left rafter/tie beam
    [2,         1,          p33],           #left rafter/king post
    [3,         0,          p06],           #right rafter/tie beam
    [3,         1,          p33],           #right rafter/king post
    [4,         1,          p13],           #left strut/king post
    [2,         4,          p22],           #left strut/left rafter
    [5,         1,          p13],           #right strut/king post
    [3,         5,          p24],           #right strut/right rafter
]

for data in connection_data :
    first_part = mdl.find_part_by_name(name='name'+str(connection_data[0]))
    second_part = mdl.find_part_by_name(name='name'+str(connection_data[1]))
    common_point = data[2]
    connector = first_part.create_connector_node(second_part.find_closest_nodes_to_point(common_point, 1, single=True))
    mdl.add_connector(RigidLinkConnector(nodes=connector, dofs="beam"))

# Set boundary conditions and insure out-of-plane stability

tie_beam = mdl.find_part_by_name(name='Part0')
king_post = mdl.find_part_by_name(name='Part1')
mdl.add_pin_bc(nodes=[tie_beam.find_closest_nodes_to_point(p00, 1, single=True)])
mdl.add_rollerXY_bc(nodes=[tie_beam.find_closest_nodes_to_point(p06, 1, single=True)])
mdl.add_rollerYZ_bc(nodes=[king_post.find_closest_nodes_to_point(p33, 1, single=True)])

# Visualize the model
# mdl.show()

#========================================
# PROBLEM
#========================================
# Define the problem
prb = mdl.add_problem(problem=Problem(name="Roof_Timber"))
# define a step
stp = prb.add_step(StaticStep())
stp.combination = LoadCombination.SLS()

# Add a node pattern to apply a load on node n2
rafter_left = mdl.find_part_by_name(name='Part2')
rafter_right = mdl.find_part_by_name(name='Part3')
stp.add_uniform_node_load(
    nodes=rafter_left.nodes + rafter_right.nodes,
    z=-1000 *units.kN,
    load_case="LL",
)

# Define the outputs
stp.add_outputs([DisplacementFieldResults, ReactionFieldResults])

#========================================
# PROBLEM
#========================================
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

#========================================
# RESULTS AND VISUALIZATION
#========================================

# # Show deformed shape
#Compas Viewer
stp.show_deformed(scale_results=10, show_original=0.3, show_loads=0.01)

# TODO Visualization of line element in VEDO
# viewer = ModelViewer(mdl)
# viewer.add_node_field_results(
#     stp.displacement_field, draw_cmap="viridis", draw_vectors=100
# )
# viewer.show()

#Reactions
react = stp.reaction_field
print("Max/Min reaction forces in X direction [kN]: ", 
      react.get_limits_component('x')[0].magnitude, "/",
      react.get_limits_component('x')[1].magnitude)
print("Max/Min reaction forces in Y direction [kN]: ", 
      react.get_limits_component('y')[0].magnitude, "/",
      react.get_limits_component('y')[1].magnitude)
print("Max/Min reaction forces in Z direction [kN]: ", 
      react.get_limits_component('z')[0].magnitude, "/",
      react.get_limits_component('z')[1].magnitude)

#TODO see the section forces