from compas_fea2.model import Model
from compas_fea2.model import DeformablePart
from compas_fea2.model import Node
from compas_fea2.model import ElasticIsotropic
from compas_fea2.model import CircularSection
from compas_fea2.model import BeamElement
from compas_fea2.model import GeneralBC

from compas_fea2.problem import Problem, StaticStep
from compas_fea2.problem import PointLoad 


### --------------------- MODEL --------------------- ###

mdl = Model()

### --------------------- PARTS --------------------- ###

prt1 = DeformablePart(name="Part1")
prt2 = DeformablePart(name="Part2")
prt3 = DeformablePart(name="Part3")

mdl.add_parts([prt1, prt2, prt3])

nodes_part1 = [Node(xyz=[x, 0.0 ,z]) for z in range(2) for x in range(0, 3, 2)]

n1p2 = Node(xyz=[1.0, 0.0, 0.0], name="Part2Node1")
n2p2 = Node(xyz=[1.0, 0.0, 1.0], name="Part2Node2")

n1p3 = Node(xyz=[1.0, 0.0, 1.0])
n2p3 = Node(xyz=[2.0, 0.0, 1.0])
n3p3 = Node(xyz=[1.5, 0.0, 2.0])

mat = ElasticIsotropic(E=210_000_000_000, v=0.2, density=7800)
sec = CircularSection(r=0.1, material=mat)

b1p1 = BeamElement(nodes=[nodes_part1[0], nodes_part1[1]], section=sec, name="Part1Element1")
b2p1 = BeamElement(nodes=[nodes_part1[1], nodes_part1[3]], section=sec, name="Part1Element2")
b3p1 = BeamElement(nodes=[nodes_part1[3], nodes_part1[2]], section=sec, name="Part1Element3")
b4p1 = BeamElement(nodes=[nodes_part1[2], nodes_part1[0]], section=sec, name="Part1Element4")
prt1.add_elements([b1p1, b2p1, b3p1, b4p1])

b1p2 = BeamElement(nodes=[n1p2, n2p2], section=sec, name="Part2Element1")
prt2.add_elements([b1p2])

b1p3 = BeamElement(nodes=[n1p3, n2p3], section=sec, name="Part1Element1")
b2p3 = BeamElement(nodes=[n2p3, n3p3], section=sec, name="Part1Element1")
b3p3 = BeamElement(nodes=[n3p3, n1p3], section=sec, name="Part1Element1")
prt3.add_elements([b1p3, b2p3, b3p3])

### ----------------------- BCs ----------------------- ###

FirstNodeBC = prt1.find_nodes_by_location(point=(0,0,0), distance=0.01)
SecondNodeBC = prt1.find_nodes_by_location(point=(2,0,0), distance=0.01)

bc_general = GeneralBC(x=True, y=True, z=True) # this is a pin
mdl.add_bcs(bc=bc_general, nodes=[SecondNodeBC[0]])

mdl.add_pin_bc(nodes=[FirstNodeBC[0]])

### ---------------------- PROBLEM ---------------------- ###


prb = Problem(name="PointLoads and Displacement")

### --- TAKE-AWAY 1: Instantiate a STATICSTEP (stp1) ------------------------------------------------------------ ###
### ---------------- Add the step to the PROBLEM ---------------------------------------------------------------- ###

stp1 = StaticStep()
prb.add_step(step=stp1)

### --- TAKE-AWAY 2: Define a vertical PointLoad to the third Node of the third Part (n3p3) --------------------- ###
### ---------------  Define an horizontal PointLoad to the fourth Node on the first Part (n4p1founded) ---------- ###

stp1.add_point_load(n3p3, z=-1000)

n4p1found = prt1.find_nodes_by_location(point=(0,0,1.0), distance=0.1)
stp1.add_point_load(n4p1found[0], x=1000)

### --- TAKE-AWAY 3: Visualize the Problem ---------------------------------------------------------------------- ###

prb.summary()
#prb.show(draw_loads=0.0005)

mdl.summary()
#mdl.show(draw_bcs=0.1)