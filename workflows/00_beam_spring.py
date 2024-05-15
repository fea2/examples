import os
from pprint import pprint 
from compas.geometry import Point, Line
from compas.colors import ColorMap, Color

import compas_fea2
from compas_fea2.model import Model, DeformablePart, SpringElement
from compas_fea2.model import ElasticIsotropic, RectangularSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_sofistik')
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

lines = [Line(Point(c*length/n, 0, z), Point((c+1)*length/n, 0, z)) for c in range(n) for z in (0, 1000)]
prt= DeformablePart.from_compas_lines(lines, section=sec)



mdl.add_part(prt)

A1 = mdl.find_nodes_around_point(Point(0, 0, 0), distance=1)
A2 = mdl.find_nodes_around_point(Point(0, 0, 1000), distance=1)
B1 = mdl.find_nodes_around_point(Point(length, 0, 0), distance=1)
B2 = mdl.find_nodes_around_point(Point(length, 0, 1000), distance=1)
prt.add_element(SpringElement(B1+B2, section=None))

mdl.add_fix_bc(A1+A2)

# mdl.summary()
# mdl.show(draw_bcs=1)

# PROBLEM
prb = mdl.add_problem(name='SLS')
stp = prb.add_static_step(nlgeom='nonl')
stp.combination = LoadCombination.SLS()
stp.add_node_pattern(nodes=B2,
                      z=1*units.kN,
                      load_case="LL")
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF'])
stp.add_output(fout)

#ANALYSIS
# mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
mdl.analyse(problems=[prb], path=TEMP, verbose=True)
