from compas_fea2.model import Model
from compas_fea2.model import DeformablePart
from compas_fea2.model import Node

from compas.geometry import Point

### --------------------- MODEL --------------------- ###

mdl = Model()

### --------------------- PARTS --------------------- ###

prt1 = DeformablePart(name="Part1")
prt2 = DeformablePart(name="Part2")
prt3 = DeformablePart(name="Part3")

mdl.add_parts([prt1, prt2, prt3])

### --- TAKE-AWAY 1: Point and Node are two different concepts. ---------------------------------------------------------------------------------------- ###
### ---------------- A POINT is a geometrical entity defined by its coordinates. ----------------------------------------------------------------------- ###
### ---------------- A NODE is an entity of the Model, it belongs to a specific PART and it is defined by its coordinates and many other attributes. --- ###
### ---------------- The same NODE can't be assign to two differents parts. The NODES n2p2 and n1p3 share the same location in the space, -------------- ### 
### ---------------- but n2p2 is assigned to the part 2 and n1p3 is assigned to the part 3 ------------------------------------------------------------- ###

nodes_part1 = [Node(xyz=[x, 0.0 ,z]) for z in range(2) for x in range(0, 3, 2)]

n1p2 = Node(xyz=[1.0, 0.0, 0.0], name="Part2Node1")
n2p2 = Node(xyz=[1.0, 0.0, 1.0], name="Part2Node2")

n1p3 = Node(xyz=[1.0, 0.0, 1.0])
n2p3 = Node(xyz=[2.0, 0.0, 1.0])
n3p3 = Node(xyz=[1.5, 0.0, 1.866]) #to change z value

### --- TAKE-AWAY 2: Nodes belongs to a specific Part, which means that after the instantiation is required -------------------------------------------- ###
### ---------------- to add them to the Part. In the next step we will see how it is not always needed ------------------------------------------------- ###

prt1.add_nodes(nodes_part1)
prt2.add_nodes([n1p2, n2p2])
prt3.add_nodes([n1p3, n2p3, n3p3])

### --- TAKE-AWAY 3: The method summary provide an overview of the Model ------------------------------------------------------------------------------- ###

mdl.summary()

