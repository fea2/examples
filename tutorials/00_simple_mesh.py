import os
from pathlib import Path
import random

from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])

mdl = Model(name='simple_frame')

lx = 10_000
ly = 10_000
nx = 10
ny = 10
mesh = Mesh.from_meshgrid(lx, nx, ly, ny)


mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = RectangularSection(w=100, h=200, material=mat)
prt = DeformablePart.frame_from_compas_mesh(mesh, sec)

mdl.add_part(prt)

fixed_nodes = [prt.find_node_by_key(vertex) for vertex in list(filter(lambda v: mesh.vertex_degree(v)==2, mesh.vertices()))]
mdl.add_fix_bc(nodes=fixed_nodes)

# DEFINE THE PROBLEM
# define a step
stp = StaticStep()
stp.combination = LoadCombination.ULS()
pt = prt.find_node_by_key(random.choice(list(filter(lambda v: mesh.vertex_degree(v)!=2, mesh.vertices()))))
stp.add_node_pattern(nodes=[pt],
                      z=-10*units.kN,
                      load_case='LL')
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF', 'S'])
stp.add_output(fout)

# set-up the problem
prb = Problem('00_simple_problem', mdl)
prb.add_step(stp)
prb.summary()

mdl.add_problem(problem=prb)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
# prb.show_displacements_contour(stp, scale_results=10, component=None, show_bcs=0.5)
prb.show_stress_contour(stp, scale_results=1e-6, component=None, show_bcs=0.5)

mdl.to_cfm(os.path.join(DATA, 'simple_frame.cfm'))