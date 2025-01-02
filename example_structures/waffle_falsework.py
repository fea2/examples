import os

from compas.datastructures import Mesh
from compas.geometry import Scale, Plane

from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput, LoadCombination

from compas_fea2.units import units
units = units(system='SI_mm')
compas_fea2.VERBOSE = False

# chage this to the backend implementation of your choice
# compas_fea2.set_backend('compas_fea2_abaqus')
# compas_fea2.set_backend('compas_fea2_sofistik')
compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
DATA = os.sep.join(HERE.split(os.sep)[:-1]+['data'])
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])

mdl = Model(name='simple_waffle')

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

settings = {       
    'target_mesh_size': 1,
    'mesh_size_at_vertices': None,
    'target_point_mesh_size': None,
    'meshsize_max': 100,
    'meshsize_min': 100,
    'rigid': False,
    "material": mat,
    "section": sec,
}
# Define a deformable part using the mesh geometry and the mechanical properties
prt = DeformablePart.from_step_file(os.path.join(DATA, 'waffle_hilo_B.stp'), **settings)
# from_gmsh(gmshModel=model, section=sec)
mdl.add_part(prt)

# fix the base
bottom_plane = Plane([0, 0, 0], [0, 0, 1])
fixed_nodes = prt.find_nodes_on_plane(bottom_plane)
mdl.add_fix_bc(nodes=fixed_nodes)
mdl.to_cfm(os.path.join(TEMP, 'waffle_hilo_B.cfm'))

mdl = Model.from_cfm(os.path.join(TEMP, 'waffle_hilo_B.cfm'))
prt = list(mdl.parts)[0]
# mdl.show(draw_bcs=0.1)

# DEFINE THE PROBLEM

# Initialize a problem
prb = mdl.add_problem(name='SLS')
# Initialize a step
stp = prb.add_static_step()
# Create a load combination
stp.combination = LoadCombination.SLS()
# Add the loads
stp.add_gravity_load_pattern([prt], g=9.81*units("m/s**2"), load_case="DL")

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'],
                   element_outputs=['S2D', 'SF'])
stp.add_output(fout)

prb.summary()
# prb.show(draw_bcs=0.1, draw_loads=0.1)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
disp = prb.displacement_field 
react = prb.reaction_field
stress = prb.stress_field

results_summary = {
    "Total volume ": mdl.volume,
    "Total weight ": mdl.volume*78/10**6*9.81/10,
    "Total vertical reaction": prb.get_total_reaction(stp)[0].z,
    "Reactions check": True if abs(round(prb.get_total_reaction(stp)[0].z-mdl.volume*23.5/10**6*9.81/10,0)) < 10 else False,
}
