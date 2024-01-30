import os
from pathlib import Path
import random

from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput

from compas_fea2.units import units
units = units(system='SI_mm')

# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')
# compas_fea2.set_backend('abaqus')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])

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
step_1 = StaticStep()
pt = prt.find_node_by_key(random.choice(list(filter(lambda v: mesh.vertex_degree(v)!=2, mesh.vertices()))))
step_1.add_node_load(nodes=[pt],
                      z=-10*units.kN)
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['SF'])
step_1.add_output(fout)
# hout = HistoryOutput('hout_test')

# set-up the problem
prb = Problem('00_simple_problem', mdl)
prb.add_step(step_1)
prb.summary()

mdl.add_problem(problem=prb)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
prb.show_deformed(scale_factor=1000, draw_loads=0.1, draw_bcs=0.25)