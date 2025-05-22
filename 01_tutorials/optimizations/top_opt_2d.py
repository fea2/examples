import os
from math import pi
from compas.datastructures import Mesh
from compas_gmsh.models import MeshModel

import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import ElasticIsotropic, ShellSection
from compas_fea2.problem import StaticStep, LoadCombination
from compas_fea2.results import StressFieldResults
import numpy as np
from compas_fea2.units import units
import matplotlib.pyplot as plt

# ==============================================================================
# UNITS: SI-mm
# ==============================================================================
units = units(system="SI_mm")
compas_fea2.set_backend("compas_fea2_abaqus")

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-2] + ["temp"])
os.makedirs(TEMP, exist_ok=True)

# ==============================================================================
# GEOMETRY & MESH
# ==============================================================================
lx = (1 * units.m).to_base_units().magnitude  # 1000 mm
ly = (30 * units.cm).to_base_units().magnitude  # 300 mm

plate = Mesh.from_polygons([[[0, 0, 0], [lx, 0, 0], [lx, ly, 0], [0, ly, 0]]])
plate = plate.rotated(pi / 2, [1, 0, 0])  # Rotate for convenience
model = MeshModel.from_mesh(plate, targetlength=30)

model.heal()
model.refine_mesh()
model.generate_mesh(2)

# ==============================================================================
# MODEL SETUP
# ==============================================================================
mdl = Model(name="beam_shell_tri")

base_E = 210000.0  # Elastic modulus [MPa] (210 GPa in N/mm²)
mat = ElasticIsotropic(E=base_E * units("MPa"), v=0.2, density=7800 * units("kg/m**3"))
sec = ShellSection(material=mat, t=5 * units.mm)

prt = Part.from_gmsh(gmshModel=model, section=sec, name="beam")
mdl.add_part(prt)

# BOUNDARY CONDITIONS
tol_x = 1e-6
for node in prt.nodes:
    if abs(node.x) < tol_x:
        mdl.add_fix_bc(nodes=[node])

# Overwrite each element's section (each gets its own material instance)
for element in prt.elements:
    mat_el = ElasticIsotropic(
        E=base_E * units("MPa"), v=0.2, density=7800 * units("kg/m**3")
    )
    sec_el = ShellSection(material=mat_el, t=5 * units.mm)
    element.section = sec_el

# ==============================================================================
# OPTIMIZATION PARAMETERS
# ==============================================================================
volfrac = 0.5  # Target volume fraction
num_elements = len(prt.elements)
densities = np.ones(num_elements) * volfrac  # Initial densities
penalty = 3.0  # Penalization factor for SIMP
move = 0.2  # Move limit
tol = 1e-4  # Convergence tolerance for densities
num_iterations = 100


# ==============================================================================
# OPTIMALITY CRITERIA FUNCTION
# ==============================================================================
def optimality_criteria_update(
    densities, strain_energy_density, volfrac, penalty, move
):
    """
    Update densities using the Optimality Criteria (OC) method without thresholds.
    """
    l1, l2 = 0.0, np.sum(
        strain_energy_density
    )  # Dynamic bounds for Lagrange multiplier
    target_vol = volfrac * len(densities)
    oc_iter = 0
    max_oc_iterations = 50

    while (l2 - l1) > tol and oc_iter < max_oc_iterations:
        lmid = 0.5 * (l1 + l2)

        # Compute updated densities
        sqrt_arg = lmid / strain_energy_density
        ideal_densities = densities * np.sqrt(sqrt_arg)

        # Apply move limit
        new_densities = np.clip(
            np.clip(ideal_densities, densities - move, densities + move), 0.0, 1.0
        )

        # Evaluate current volume
        current_vol = np.sum(new_densities)
        if current_vol > target_vol:
            l1 = lmid  # Too much volume, increase penalty
        else:
            l2 = lmid  # Too little volume, decrease penalty

        oc_iter += 1

    change = np.max(np.abs(new_densities - densities))
    return new_densities, change


# ==============================================================================
# OPTIMIZATION LOOP
# ==============================================================================
for i in range(num_iterations):
    print("\n==========================")
    print(f"Iteration {i+1}/{num_iterations}")
    print("==========================")

    # Create new problem & step
    prb = mdl.add_problem(name=f"iter-{i}")
    stp = StaticStep(
        name=f"iter-{i}", system="SparseGeneral", test="NormDispIncr 1.0e-1 10"
    )
    prb.add_step(stp)
    stp.combination = LoadCombination.SLS()

    # Apply load to the right edge
    loaded_nodes = [n for n in prt.nodes if abs(n.x - lx) < tol_x]
    load_value = -(1.0 / len(loaded_nodes)) * units.kN
    stp.add_uniform_node_load(nodes=loaded_nodes, z=load_value, load_case="LL")
    stp.add_output(StressFieldResults)

    # Update material properties based on current densities
    for idx, element in enumerate(prt.elements):
        element.section.material.E = densities[idx] ** penalty * base_E

    # Perform FEA analysis
    prb.analyse_and_extract(path=TEMP,  output=True)

    # Extract & normalize strain energy density
    strain_energy_density = np.zeros(num_elements)
    for idx, element in enumerate(prt.elements):
        try:
            results = element.stress_results(stp)
            mid_stress = results.mid_plane_stress_result
            sed_val = max(
                mid_stress.strain_energy_density, 1e-12
            )  # Avoid division by zero
            strain_energy_density[idx] = sed_val
        except (AttributeError, ValueError) as e:
            print(f"Error extracting stress results for element {idx}: {e}")
            strain_energy_density[idx] = 1e-12

    # Normalize SED to avoid unit-based inconsistencies
    strain_energy_density /= np.max(strain_energy_density)

    # Call the OC update
    new_densities, change = optimality_criteria_update(
        densities, strain_energy_density, volfrac, penalty, move
    )
    densities = new_densities

    # Print iteration stats
    print(f"Iteration {i+1}: Change in densities = {change:.4e}")
    print(f"Iteration {i+1}: Volume fraction = {np.sum(densities) / num_elements:.3f}")

    if change < tol:
        print("Optimization converged.")
        break

# ==============================================================================
# FINAL DENSITY PLOT
# ==============================================================================
plt.figure(figsize=(10, 6))
plt.plot(range(len(densities)), densities, "bo-", label="Final Densities")
plt.xlabel("Element Index")
plt.ylabel("Density")
plt.title("Final Optimized Densities")
plt.legend()
plt.show()
