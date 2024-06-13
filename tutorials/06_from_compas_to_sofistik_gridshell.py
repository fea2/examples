import os
from pathlib import Path
import random

import compas
from compas.datastructures import Mesh

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_sofistik')

my_cwd = os.getcwd()
DIR = os.path.join(my_cwd, "data", "tutorials", "06_from_compas_to_sofistik_gridshell")

mdl = Model(name='my_model')

lx = (10*units.m).to_base_units().magnitude
ly = (10*units.m).to_base_units().magnitude
nx = 10
ny = 10
mesh = Mesh.from_meshgrid(lx, nx, ly, ny)

mat = ElasticIsotropic(E=(210*units.GPa).to_base_units().magnitude, 
                       v=0.2, 
                       density=(7800*units("kg/m**3")).to_base_units().magnitude)
sec = RectangularSection(w=100, h=200, material=mat)
prt = DeformablePart.frame_from_compas_mesh(mesh, sec)

mdl.add_part(prt)

fixed_nodes = [prt.find_node_by_key(vertex) for vertex in list(filter(lambda v: mesh.vertex_degree(v)==2, mesh.vertices()))]
mdl.add_fix_bc(nodes=fixed_nodes)

step_1 = StaticStep()
pt = prt.find_node_by_key(random.choice(list(filter(lambda v: mesh.vertex_degree(v)!=2, mesh.vertices()))))
step_1.add_point_load(nodes=[pt],
                      z=-(10*units.kN).to_base_units().magnitude)

#The name "from_compas_to_sofistik_gridshell" will be the .dat file name to run in sofistik
prb = Problem('from_compas_to_sofistik_gridshell', mdl)
prb.add_step(step_1)
prb.summary()

mdl.add_problem(problem=prb)

#the location where the file from_compas_to_sofistik_gridshell.dat will be saved.
prb.path = DIR

#command to generate the .dat file.
prb.write_input_file()