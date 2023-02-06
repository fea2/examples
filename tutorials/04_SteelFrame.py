import compas_fea2

from compas_fea2.model import Model 
from compas_fea2.model import DeformablePart 

from compas_fea2.model import Node 
from compas_fea2.model import BeamElement 
from compas_fea2.model import ISection
from compas_fea2.model import ElasticIsotropic 

from compas_fea2.problem import Problem, StaticStep
from compas_fea2.problem import PointLoad 

from compas_fea2.model import PinnedBC


mdl = Model()
prt = DeformablePart()

n1 = Node(xyz=[0.0, 0.0, 0.0], name="FirstNode")
n2 = Node(xyz=[0.0, 0.0, 6.0])
n3 = Node(xyz=[7.6, 0.0, 6.0])
n4 = Node(xyz=[7.6, 0.0, 0.0])

mat_steel = ElasticIsotropic(E=210_000_000_000, v=0.2, density=7800)
sec_column = ISection(w=0.26, h=0.26, tw=0.0075, tf=0.0125, material=mat_steel)              
sec_beam = ISection(w=0.2, h=0.5, tw=0.012, tf=0.016, material=mat_steel)                     

column1 = BeamElement(nodes=[n1, n2], section=sec_column)
column2 = BeamElement(nodes=[n4, n3], section=sec_column)
beam = BeamElement(nodes=[n2, n3], section=sec_column)
prt.add_elements([column1, column2, beam])

mdl.add_parts([prt])
mdl.add_bcs(bc=PinnedBC(), nodes=[n1, n4])

#old
#prb = Problem(model=mdl)

#trying new
prb = Problem()
#prb.model = mdl

stp = StaticStep()
stp.add_point_load(n2, x=1000)

prb.add_step(step=stp)
#mdl.show()
#mdl.show(draw_loads=0.001)

mdl.summary()
