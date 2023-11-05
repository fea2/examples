import os

from compas.datastructures import Mesh
from compas.geometry import Scale, Plane

from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
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

mdl = Model(name='simple_waffle')

# Get the waffle geometry
mesh = Mesh.from_obj(os.path.join(DATA, 'simple_waffle', 'simple_waffle.obj'))
S = Scale.from_factors([1000.] * 3)
mesh.transform(S)

# Use GMSH to discretize the geometry
print('generating discretization...')
model = MeshModel.from_mesh(mesh, targetlength=500)
model.heal()
model.generate_mesh(2)
compas_mesh = model.mesh_to_compas()
print("discretization complete!")

# Define mechanical properties and panel thickness
E = 10*units.GPa
v = 0.2
rho = 500*units("kg/m**3")
t = 20*units.mm

# Define material and section (here it is the same for each element)
mat = ElasticIsotropic(E=E.to_base_units().magnitude,
                       v=v,
                       density=rho.to_base_units().magnitude)
sec = ShellSection(t=t.to_base_units().magnitude,
                   material=mat)

# Define a deformable part using the mesh geometry and the mechanical properties
prt = DeformablePart.shell_from_compas_mesh(mesh=compas_mesh, section=sec)
# from_gmsh(gmshModel=model, section=sec)
mdl.add_part(prt)

# fix the base
bottom_plane = Plane([0, 0, 0], [0, 0, 1])
fixed_nodes = prt.find_nodes_on_plane(bottom_plane)
mdl.add_fix_bc(nodes=fixed_nodes)
# mdl.show()

# DEFINE THE PROBLEM
prb = Problem('point_loads', mdl)
# define a Linear Elastic Static analysis
step_1 = StaticStep()
# Define the loads
top_plane = Plane([0, 0, 1000], [0, 0, 1])
loaded_nodes = prt.find_nodes_on_plane(top_plane)
step_1.add_point_load(nodes=loaded_nodes,
                      z=-(10*units.N).to_base_units().magnitude)
# step_1.add_gravity_load(g=9810)
# decide what information to save
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S'])
step_1.add_output(fout)

# Add the step to the problem
prb.add_step(step_1)
prb.summary()

# Add the problem to the model
mdl.add_problem(problem=prb)
mdl.summary()

# Run the analysis
# mdl.analyse(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
print(f'Analysis results saved in {prb.path}')
