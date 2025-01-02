import os
from pprint import pprint 
from compas.geometry import Point, Line
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node, BeamElement
from compas_fea2.model import ElasticIsotropic, RectangularSection, ISection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')
compas_fea2.VERBOSE = False

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2]+['temp'])


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

sec = ISection(w=25*units.cm, h=250*units.mm, tw=20*units.mm, tf=4*units.cm, material=mat)

length = 1000
n = 10

lines = [Line(Point(c*length/n, 0, 0), Point((c+1)*length/n, 0, 0)) for c in range(n)]
prt= DeformablePart.from_compas_lines(lines, section=sec)
mdl.add_part(prt)

A = mdl.find_nodes_around_point(Point(0, 0, 0), distance=1)[0]
B = mdl.find_nodes_around_point(Point(length, 0, 0), distance=1)[0]

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
prb.path=TEMP
prb.write_input_file(path=TEMP)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
prb.show_displacements_contour()


disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field

print(f"Total volume : {mdl.volume*units('mm**3').to('m**3'):P}")

print(prb.results_db.fields)
print("\nDisplacement vector with min component along 2:")
print(disp.get_min_component(2, step=stp).vector)
print("Displacement vector at B:")
print(disp.get_results_at_point(B.point, 1, steps=[stp])[stp][0].vector)

print("\nDisplacement vectors with max and min component along 2:")
for r in disp.get_limits_component(2, step=stp):
    print(r.vector)

print("\nAbsolute max and min displacement vectors :")
for r in disp.get_limits_absolute(step=stp):
    print(r.vector)
