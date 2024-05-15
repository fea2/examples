import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import CircularSection, ElasticIsotropic, BeamElement
from compas_fea2.units import units
units = units(system='SI_mm')

from compas.geometry import Frame


mdl = Model(name='my_model')
prt = DeformablePart(name='my_part')
n1 = Node(xyz=[0,0,0])
n2 = Node(xyz=[1,0,0])

mat = ElasticIsotropic(E=210*units("GPa"), v=0.2, density=7800*units("kg/m**3"))
sec = CircularSection(r=10*units.cm, material=mat)
beam = BeamElement(nodes=[n1,n2], section=sec, frame=Frame.from_points(n1.point, n2.point, [0, 1, 0]))

prt.add_element(beam)
print(prt)

mdl.add_part(part=prt)

mdl.summary()
mdl.show()