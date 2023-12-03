from compas_fea2.units import units
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import ElasticIsotropic
from compas_fea2.model import RectangularSection, CircularSection, PipeSection, ISection, BoxSection, AngleSection
from compas_fea2.model import BeamElement
from compas_fea2.model import GeneralBC

from compas_fea2.problem import Problem
from compas_fea2.problem import StaticStep
from compas_fea2.problem import PointLoad, LineLoad

import os

### --- TAKE-AWAY: How to set the backend sofistik ------------------------------------------------------------ ###
### -------------- How to name the .dat input file ------------------------------------------------------------ ###
### -------------- How to define location for the generation of the .dat input file --------------------------- ###
### -------------- Use the 'write_inpute_file()' command ------------------------------------------------------ ###

import compas_fea2
compas_fea2.set_backend('compas_fea2_sofistik')

my_cwd = os.getcwd()
DIR = os.path.join(my_cwd, "data", "tutorials", "05_from_compas_to_sofistik_simple_beam")

units = units(system='SI_mm')

mdl = Model()
prt = DeformablePart()
mdl.add_part(part=prt)

nodes = [Node(xyz=[0, i*1000, 0]) for i in range(7)]
prt.add_nodes(nodes=nodes)

mat = ElasticIsotropic(E=(210*units.GPa).to_base_units().magnitude,
                       v=0.2,
                       density=(7800*units("kg/m**3")).to_base_units().magnitude)
sec = RectangularSection(w=100, h=200, material=mat)

beams_elements = [BeamElement(nodes=[nodes[i], nodes[i+1]], section=sec) for i in range(len(nodes) - 1)]
prt.add_elements(beams_elements)

mdl.add_fix_bc(nodes=[nodes[0], nodes[len(nodes)-1]])

#The name "from_compas_to_sofistik_simple_beam" will be the .dat file name to run in sofistik
prb = Problem("from_compas_to_sofistik_simple_beam", mdl)
stp = StaticStep()
prb.add_step(step=stp)

NodesPointLoad2 = [prt.find_nodes_by_location(point=(0, 1000, 0), distance=0.01)[0]]
NodesPointLoad3 = [prt.find_nodes_by_location(point=(0, 2000, 0), distance=0.01)[0]]
NodesPointLoad4 = [prt.find_nodes_by_location(point=(0, 3000, 0), distance=0.01)[0]]

stp.add_point_load(nodes=NodesPointLoad2, z=5.0)
stp.add_point_load(nodes=NodesPointLoad4, x=2.0)

mdl.add_problem(problem=prb)

#the location where the file from_compas_to_sofistik_simple_beam.dat will be saved.
prb.path = DIR

#command to generate the .dat file.
prb.write_input_file()
