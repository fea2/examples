from compas_fea2.model import ISection, Steel

isection = ISection.IPE160(material=Steel.S355())
ishape = isection.shape.plot()

isection = ISection.HEA200(material=Steel.S355())
ishape = isection.shape.plot()

isection = ISection.HEB650(material=Steel.S355())
ishape = isection.shape.plot()
