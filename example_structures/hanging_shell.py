import os
from pathlib import Path
import random

import compas
from compas.datastructures import Mesh
from compas.geometry import Translation, Vector, Scale, Frame, Plane

from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput

from compas_fea2.units import units
units = units(system='SI_mm')

# chage this to the backend implementation of your choice
# compas_fea2.set_backend('compas_fea2_abaqus') 
# compas_fea2.set_backend('compas_fea2_sofistik') 
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])

mdl = Model(name='hanging_shell')

# Get the geometry from the obj file (and scale it to mm)
mesh = Mesh.from_obj(os.path.join(DATA, 'hanging_shell', 'tofea_f.obj'))
S = Scale.from_factors([1000.] * 3)
mesh.transform(S)

# Use GMSH to discretize the geometry
print('generating discretization...')
model = MeshModel.from_mesh(mesh, targetlength=100)
model.heal()
model.generate_mesh(2)
compas_mesh = model.mesh_to_compas()
print("discretization complete!")

# Define mechanical properties and panel thickness
E = 10*units.GPa
v = 0.2
rho = 1500*units("kg/m**3")
t = 20*units.mm

# Define material and section (here it is the same for each element)
mat = ElasticIsotropic(E=E, 
                       v=v, 
                       density=rho)
sec = ShellSection(t=t, 
                   material=mat)

# Define a deformable part using the mesh geometry and the mechanical properties
prt = DeformablePart.shell_from_compas_mesh(mesh=compas_mesh, section=sec)
mdl.add_part(prt)

# fix the base
bottom_plane = Plane([0,0,0], [0,0,1])
fixed_nodes = prt.find_nodes_on_plane(bottom_plane)
mdl.add_pin_bc(nodes=fixed_nodes)
# mdl.show()

# DEFINE THE PROBLEM
prb = Problem('gravity', mdl)

# define a Linear Elastic Static analysis
step_1 = StaticStep()

# Define the loads
step_1.add_gravity_load(g=9.81*units("m/s**2"), z=-1.)
# step_1.add_point_load(prt.nodes, z=-(10*units.kN).to_base_units().magnitude)

# decide what information to save
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D', 'SF'])
step_1.add_output(fout)

# Add the step to the problem
prb.add_step(step_1)
# prb.summary()

# Add the problem to the model
mdl.add_problem(problem=prb)
# mdl.summary()

# ANALYSIS and RESULTS
# Run the analysis
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
print(f'Analysis results saved in {prb.path}')
# prb.show_nodes_field_contour('S', 3)
# prb.show_nodes_field_vector('U', vector_sf=500)
# prb.show_elements_field_vector('S', vector_sf=50)
prb.show_elements_field_vector('S2D', vector_sf=50, draw_bcs=0.5, draw_loads=0.1)
# prb.show_deformed(scale_factor=500, draw_bcs=0.5, draw_loads=1000, opacity=1., original=0.25)
# prb.show_deformed()