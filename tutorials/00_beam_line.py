import os

from compas.datastructures import Mesh
from compas.colors import ColorMap, Color
from compas.geometry import Line

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, RectangularSection, ISection, BeamSection
from compas_fea2.model.shapes import Rectangle
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination
from compas_fea2.results import DisplacementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


l = (0.1*units.m).to_base_units().magnitude
n = 50

lines = [Line([i*l, 0, 0], [(i+1)*l, 0, 0]) for i in range(n-1)]

# Initialize model
mdl = Model(name='beam_line')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ISection(w=100, h=190, tw=2, tf=2, material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt= DeformablePart.from_compas_lines(lines, section=sec, name='beam')
mdl.add_part(prt)

# Set boundary conditions
mdl.add_fix_bc(nodes=prt.find_closest_nodes_to_point([0, 0, 0], distance=0.1))

mdl.summary()
# mdl.show(show_bcs=0.2)

prb = mdl.add_problem(name='beam_line_Fz')
stp = prb.add_static_step()
stp.combination = LoadCombination.SLS()

# Add the load
stp.add_node_pattern(nodes=prt.find_closest_nodes_to_point([(n-1)*l, 0, 0], distance=0.1), x=-1*units.kN, z=-1*units.kN, load_case='LL')

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF'])
stp.add_output(fout)

# prb.show(show_loads=0.1, show_bcs=0.1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=os.path.join(TEMP, prb.name), verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
    
print("Displacement vector with min z component: ", disp.get_min_result(3, stp).vector)
print("Min z component of the displacement field [mm]: ", disp.get_min_component(3, stp))

prb.show_deformed(scale_results=10, show_original=0.1, show_bcs=0.1)