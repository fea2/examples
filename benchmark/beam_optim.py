import os

from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, DeformablePart
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput
from compas_fea2.results import NodeFieldResults, ElementFieldResults

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')
compas_fea2.VERBOSE = False

from scipy.optimize import minimize
# from pyOpt import Optimization, ALPSO

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


def perform_fea(thicknesses, elements):
    for c, e in enumerate(elements):
        e.section = ShellSection(material=mat, t=e.section.t+thicknesses[c]) 
    mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
    stress_field = ElementFieldResults('S2D', stp)
    return [r.von_mises_stress for r in stress_field.results]
    
# Define the objective function to minimize weight while considering stress constraints
def objective_function(thicknesses, elements, current_iteration):
    # Increment the current iteration counter
    print(current_iteration)
    # Perform FEA using your library to get von Mises stresses
    vms_results = perform_fea(thicknesses, elements)
    
    # Calculate the total weight based on the thicknesses
    weight = sum([prt.volume for prt in mdl.parts])# Implement this function

    # Penalize violations of the stress constraint (e.g., using a penalty factor)
    stress_limit = 10  # Set your von Mises stress limit here
    stress_penalty = 10000  # Adjust the penalty factor as needed
    stress_violations = max(0, max(vms_results) - stress_limit)

    # Combine the weight and stress constraint in the objective
    total_objective = weight + stress_penalty * stress_violations
    current_iteration += 1
    return total_objective

# Define your custom callback function
def callback(xk):
    # This function will be called at each iteration with the current design variables (xk)
    print("Iteration {}: Thicknesses: {}".format(callback.iteration, xk))
    callback.iteration += 1

# ==============================================================================
# Make a plate mesh
# ==============================================================================
lx = (1*units.m).to_base_units().magnitude
ly = (30*units.cm).to_base_units().magnitude

plate = Mesh.from_polygons([[[0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0]]])
model = MeshModel.from_mesh(plate, targetlength=50)

model.heal()
model.refine_mesh()
model.generate_mesh(2)
# ==============================================================================
# COMPAS_FEA2
# ==============================================================================

# Initialize model
mdl = Model(name='plate')
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                       v=0.2, 
                       density=7800*units("kg/m**3"))
sec = ShellSection(material=mat, t=30*units.mm)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec, name='beam', ndm=3, implementation='shelldkgt')
prt._discretized_boundary_mesh = model.mesh_to_compas()
prt._boundary_mesh = plate
prt.bounding_box
mdl.add_part(prt)

# Set boundary conditions in the corners
for node in prt.nodes:
    if node.x == 0:
        mdl.add_fix_bc(nodes=[node])

mdl.summary()
# mdl.show(draw_bcs=0.1)

# Initialize a step
stp = StaticStep()

# Add the load
loaded_nodes = list(filter(lambda n: n.x == lx, prt.nodes))
stp.add_node_load(nodes=loaded_nodes,
                  y=-(2/len(loaded_nodes))*units.kN)

# Ask for field outputs
fout = FieldOutput(element_outputs=['S2D'])
stp.add_output(fout)

# Set-up the problem
prb = Problem('beam_shell_opti')
prb.add_step(stp)
mdl.add_problem(problem=prb)

# Analyze and extracte results to SQLite database
mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)

stress_field = ElementFieldResults('S2D', stp)
# Define initial thicknesses and bounds
initial_thicknesses = [r.element.section.t for r in stress_field.results]
elements = [r.element for r in stress_field.results]


from scipy.optimize import minimize
import numpy as np

# Define the number of elements (you need to set this)
numElements = 10


# Define the lower and upper bounds for thickness for each element
lower_bounds = [1] * numElements  # Replace with your minimum thickness
upper_bounds = [200] * numElements  # Replace with your maximum thickness

# Define the objective and constraint function
def objective_andConstraint(thickness):
    for e, t in zip(elements, thickness):
        e.section = ShellSection(t=t, material=mat)
    
    # mdl.show()
    mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
    weight = sum([p.weight for p in mdl.parts])
    print(weight)
    
    stress_field = ElementFieldResults('S2D', stp) 
    vms = [r.von_mises_stress for r in stress_field.results]
    max_stresses = [50]*len(vms)
    stress_violations = [max(0, stress - max_stress) for stress, max_stress in zip(vms, max_stresses)]
    
    return weight, stress_violations

# Define the optimization problem using scipy.optimize
result = minimize(lambda x: objective_andConstraint(x)[0], initial_thicknesses, options={'maxiter': 1, 'maxfun':1, 'disp':True})
# result = minimize(lambda x: objective_andConstraint(x)[0], initial_thicknesses, bounds=list(zip(lower_bounds, upper_bounds)), constraints={'type': 'ineq', 'fun': lambda x: objective_andConstraint(x)[1]})

# Access the optimal thickness values for each element
optimal_thickness = result.x



# # Define the objective function
# def objective(thickness):
#     # Update the existing 'mdl' instance with the new thickness values
#     for e, t in zip(elements, thickness):
#         e.section.t = t

#     # Analyze and extracte results to SQLite database
#     mdl.analyse_and_extract(problems=[prb], path=TEMP, verbose=True)
#     weight = sum([p.volume for p in mdl.parts])

#     stress_field = ElementFieldResults('S2D', stp) 
    
#     vms = [r.von_mises_stress for r in stress_field.results]
       
#     # Set the maximum stress threshold for each element (replace with your desired values)
#     max_stresses = [50]*len(vms)
    
#     # Calculate the stress constraint violations for each element
#     stress_violations = [max(0, stress - max_stress) for stress, max_stress in zip(vms, max_stresses)]
    
#     # Return weight as the objective and stress constraint violations as constraints
#     return weight, stress_violations


# # Initialize PyOpt problem
# opt_problem = Optimization('FE Analysis Optimization', objective)
# opt_problem.addObj('Weight')


# # Add stress constraint for each element
# for i in range(elements):
#     opt_problem.addCon(f'Stress_{i}', type='ineq', lower=0.0)


# min_thk = 1
# max_thk = 200
# # Add variables and bounds for each element's thickness
# for i in range(elements):
#     opt_problem.addVar(f'thickness_{i}', type='c', value=initial_thicknesses[i], 
#                        lower=min_thk, upper=max_thk)

# # Initialize the optimizer (e.g., ALPSO)
# optimizer = ALPSO()
# optimizer.setOption('SwarmSize', 20)  # Adjust as needed
# optimizer.setOption('maxOuterIter', 100)  # Adjust as needed

# # Solve the optimization problem
# optimizer(opt_problem)

# # Access the optimal thickness values for each element
# optimal_thickness = [opt_problem.solution(i) for i in range(elements)]









mdl.show(draw_bcs=0.05, drawLoads=0.1)
# # Show Results
# prb.show_nodes_field_contour('U', '2')
# prb.show_nodes_field_vector('U', vector_sf=500)
# prb.showElements_field_vector('S2D', vector_sf=1, draw_bcs=0.05, drawLoads=0.1)