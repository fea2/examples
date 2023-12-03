# COMPAS_FEA2_POEM

# A compas_fea2 Model consists of several Parts.
# A Part consists of Nodes and Elements.
# To define an Element is needed to define the Nodes and a Section.
# To define a Section is needed to define a Material.
# BCs are assigned to the Nodes
# A Problem consists of several Steps
# Exist different types of Steps
# From the Results you can get Outputs
# From the Results you can get Results

# The goal of the following examples is to provide an understanding of the main concepts in compas_fea2

from compas_fea2.model import Model

### --------------------- MODEL --------------------- ###

mdl = Model()