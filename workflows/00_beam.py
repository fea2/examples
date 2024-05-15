import os
from pprint import pprint 
from compas.geometry import Point, Line
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node, BeamElement
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')
compas_fea2.VERBOSE = False

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a beam with line elements
# ==============================================================================


# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='beam_line')
# mdl._ndf = 3
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7.8*units("ton/m**3"))

sec = RectangularSection(w=50*units.mm, h=100*units.mm, material=mat)

# Convert the gmsh model in a compas_fea2 Part
length = 5000
n = 10

lines = [Line(Point(c*length/n, 0, 0), Point((c+1)*length/n, 0, 0)) for c in range(n)]
prt= DeformablePart.from_compas_lines(lines, section=sec)
mdl.add_part(prt)

A = mdl.find_nodes_around_point(Point(0, 0, 0), distance=1)
B = mdl.find_nodes_around_point(Point(length, 0, 0), distance=1)

mdl.add_fix_bc(A)

# mdl.summary()

# PROBLEM
prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()
stp.add_node_pattern(nodes=B,
                      y=-1*units.kN,
                      load_case="LL")
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF'])
stp.add_output(fout)

#ANALYSIS
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field

# results_summary = {
#     "Total volume : {mdl.volume*units('mm**3').to('m**3')}",
#     "Total weight ": (mdl.volume*78/10**6*9.81/10/1000)*units('kN'),
#     "Total vertical reaction": prb.get_total_reaction(stp)[0].y/1000*units('kN'),
# }
print(f"Total volume : {mdl.volume*units('mm**3').to('m**3'):P}")

# print(prb.results_db.fields)
# print(disp.get_results(members=pt, steps=stp))
# print(disp.get_min_component(3, step=stp).vector)
# print(disp.get_limits_component(3, step=stp))
# print(disp.get_limits_absolute(step=stp))
# # print(disp.all_results[0])
# # print(disp.get_value_at_node(pt[0], stp))
# # print(disp.max.magnitude)
# # print(disp.min.magnitude)

# cmap = ColorMap.from_color(Color.red(), rangetype='light')
cmap = ColorMap.from_mpl('viridis')

# Show Results
prb.show_nodes_field_contour(disp, component=3, draw_reactions=0.1, draw_loads=0.1, draw_bcs=0.1, cmap=cmap)
prb.show_nodes_field_vector(disp, component=3, scale_factor=1000, draw_bcs=0.1,  draw_loads=0.1)
prb.show_deformed(scale_factor=1000, draw_bcs=0.1, draw_loads=0.1)
prb.show_stress_contours(stress_type="von_mises_stress", side="top", draw_reactions=0.1, draw_loads=0.1, draw_bcs=0.1, cmap=cmap, bounds=[0, 0.5])
prb.show_elements_field_vector(stress, vector_sf=10, draw_bcs=0.1)
