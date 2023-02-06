import compas_fea2

from compas_fea2.model import Model
from compas_fea2.model import DeformablePart
from compas_fea2.model import Node

from compas.geometry import Point

### --- TAKE-AWAY 1: understand that defining the Nodes without adding them to the DeformablePart is not enough to register them to the Model ----------------- ###

mdl = Model()
prt = DeformablePart()

mdl.add_parts([prt])

n1 = Node(xyz=[0.0, 0.0, 0.0], name="node_n1")
n2 = Node(xyz=[1.0, 0.0, 0.0], name="node_n2")
n3 = Node(xyz=[2.0, 0.0, 0.0], name="node_n3")

#Note: Since only the nodes n1 and n3 will be added to the DeformablePart, when we will use the find_nodes_by_location with the distance equal to 10, 
#it will find only n1 and n3. It means that even if n2 was instantiated as Node, it doesn't belong to the DeformablePart 'prt'. 

prt.add_nodes([n1, n3])

node_found = prt.find_nodes_by_location([1.0, 1.0, 1.0], distance=10.0)
print(*[node.name for node in node_found])




