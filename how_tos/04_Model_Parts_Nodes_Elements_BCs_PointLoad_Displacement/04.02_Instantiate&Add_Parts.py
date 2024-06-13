# After we instantiate a Model and one or more Parts, it is needed to add the Parts to the Model.
# In this example, our Model has three parts. The Model will be instantiated as mdl, while the three parts as prt1, prt2, prt3. 


### --- TAKE-AWAY 1: Instantiate few Parts and add them to the Model -------------------------- ###

from compas_fea2.model import Model
from compas_fea2.model import DeformablePart
from compas_fea2.model import Node

### --------------------- MODEL --------------------- ###

mdl = Model()

### --------------------- PARTS --------------------- ###

prt1 = DeformablePart(name="Part1")
prt2 = DeformablePart(name="Part2")
prt3 = DeformablePart(name="Part3")

mdl.add_parts([prt1, prt2, prt3])

