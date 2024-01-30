import time
import os
from random import choice
from pathlib import Path

from compas.geometry import Point, Polyline, Bezier, Vector
from compas.datastructures import Mesh, mesh_thicken
from compas.colors import ColorMap, Color
from compas.utilities import geometric_key_xy

from compas_gmsh.models import MeshModel

from compas_view2.app import App
from compas_view2.shapes import Arrow

import compas_fea2
from compas_fea2.model import Model, DeformablePart, Node
from compas_fea2.model import RectangularSection, ElasticIsotropic, SolidSection
from compas_fea2.problem import Problem, StaticStep, FieldOutput

from compas_fea2.units import units
units = units(system='SI_mm')

compas_fea2.set_backend('compas_fea2_opensees')

HERE = os.path.dirname(__file__)
TEMP = os.sep.join(HERE.split(os.sep)[:-1]+['temp'])


# ==============================================================================
# Make a plate mesh and discretize it with GMSH
# ==============================================================================
lx = (10*units.m).to_base_units().magnitude
ly = (10*units.m).to_base_units().magnitude
nx = 5
ny = 5
mesh = Mesh.from_meshgrid(lx, nx, ly, ny)
plate = mesh_thicken(mesh, 100)

poa = choice(list(set(mesh.vertices()) - set(mesh.vertices_on_boundary())))
poa_coordinates = mesh.vertex_coordinates(poa)

model = MeshModel.from_mesh(plate, targetlength=500)

# refine mesh at the point of application of the load
model.mesh_targetlength_at_vertex(poa, 50)

# refine mesh at the supports
for vertex in mesh.vertices_where({'vertex_degree': 2}):
    a = geometric_key_xy(mesh.vertex_coordinates(vertex))
    for vertex in plate.vertices():
        b = geometric_key_xy(plate.vertex_coordinates(vertex))
        if a == b:
            model.mesh_targetlength_at_vertex(vertex, 50)

model.heal()
model.generate_mesh(3)
# model.optimize_mesh(niter=10)
# model.recombine_mesh()

solid_mesh = model.mesh_to_compas()

# ==============================================================================
# COMPAS_FEA2
# ==============================================================================
# Initialize model
mdl = Model()
# Define some properties
mat = ElasticIsotropic(E=210*units.GPa, 
                    v=0.2, 
                    density=7800*units("kg/m**3"))
sec = SolidSection(material=mat)

# Convert the gmsh model in a compas_fea2 Part
prt = DeformablePart.from_gmsh(gmshModel=model, section=sec)
prt.ndf=3 # this is needed for the opensees FourNodeTetrahedron model
prt._discretized_boundary_mesh = solid_mesh
mdl.add_part(prt)

# Set boundary conditions in the corners
for vertex in mesh.vertices_where({'vertex_degree': 2}):
    location = mesh.vertex_coordinates(vertex)
    mdl.add_fix_bc(nodes=prt.find_nodes_by_location(location, distance=150))

# Initialize a step
stp = StaticStep()

# Add the load
pt = prt.find_closest_nodes_to_point(poa_coordinates, distance=150)
stp.add_node_load(nodes=pt,
                    z=-10*units.kN)

# Ask for field outputs
fout = FieldOutput(node_outputs=['U', 'RF'])
stp.add_output(fout)

# Set-up the problem
prb = Problem('01_mesh_refine')
prb.add_step(stp)
mdl.add_problem(problem=prb)

# ==============================================================================
# Viz
# ==============================================================================

viewer = App(viewmode="shaded", enable_sidebar=True, width=1600, height=900)

viewer.view.camera.rz = 0
viewer.view.camera.rx = -55
viewer.view.camera.tx = -5
viewer.view.camera.ty = -2
viewer.view.camera.distance = 10

viewer.view.camera.target = [5000, 5000, 100]
viewer.view.camera.position = [5000, 0, 5000]
viewer.view.camera.near = 10
viewer.view.camera.far = 100000
viewer.view.camera.scale = 1000
viewer.view.grid.cell_size = 1000

obj = viewer.add(solid_mesh)

poa = Point(* plate.vertex_coordinates(poa))
start = poa + Vector(0, 0, 1000)
vector = Vector(0, 0, -1000)
load = Arrow(start, vector, body_width=0.02)

viewer.add(poa, pointsize=5)
viewer.add(load, facecolor=(1, 0, 0))

@viewer.button(text="Compute")
def click():

    viewer.status(f"Analysis running...")

    # Function to be executed in the background thread.
    def compute(self, step, interval):
        self.signals.progress.emit(1)

    def on_progress(value):
        # This function will be triggered under main thread when `signals.progress.emit()` sends out value from background threads.

        # Analyze and extracte results to SQLite database
        mdl.analyse(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)
        # mdl.analyse_and_extract(problems=[prb], path=Path(TEMP).joinpath(prb.name), verbose=True)

        # Serialize
        prb.convert_results_to_sqlite(Path(TEMP).joinpath(prb.name, prb.name))
        
        
        c_symb = 'U3'
        c_index = 2

        displacements, _ = prb.get_displacements_sql()
        max_disp = prb.get_max_displacement_sql(component=c_symb)
        min_disp = prb.get_min_displacement_sql(component=c_symb)

        parts_gkey_vertex={}
        parts_mesh={}
        
        for part in prb.model.parts:
            if (mesh:= part.discretized_boundary_mesh):
                colored_mesh = mesh.copy()
                parts_gkey_vertex[part.name] = colored_mesh.gkey_key(compas_fea2.PRECISION)
                parts_mesh[part.name] = colored_mesh

        pts = []
        vectors = []
        colors = []

        for displacement in displacements:
            part = displacement['part']
            node = displacement['node']
            vector = displacement['vector']
            
            cmap = ColorMap.from_palette('hawaii') #ColorMap.from_color(Color.red(), rangetype='light') #ColorMap.from_mpl('viridis')
            color = cmap(vector[c_index], minval=min_disp['vector'][c_index], maxval=max_disp['vector'][c_index])

            pts.append(node.point)
            vectors.append(vector)
            colors.append(color)

            if part.discretized_boundary_mesh:
                if node.gkey in parts_gkey_vertex[part.name]:
                    parts_mesh[part.name].vertex_attribute(parts_gkey_vertex[part.name][node.gkey], 'color', color)

        viewer.remove(obj)
        for part in prb.model.parts:
            if part.discretized_boundary_mesh:
                viewer.add(parts_mesh[part.name], use_vertex_color=True)

        
        viewer.view.update()

    def on_result(result):
        # This function will be triggered once the background thread finishes.
        viewer.status(f"Analysis complete")
        viewer.info("Finished!")

    # `include_self=True` is provide in order to give "compute" function access to singal emit functions.
    viewer.threading(compute, args=[0.02, 0.2], on_progress=on_progress, on_result=on_result, include_self=True)


viewer.show()