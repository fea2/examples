# step_file="/Users/francesco/code/VAULTED/RFS_fea/00_stp_test_models/rfs-d.stp",  # noqa: E501
import os
import compas_fea2
from compas_fea2.model import Model, Part
from compas_fea2.model import SolidSection, Steel
from compas_fea2.units import units
from compas.datastructures import Mesh
from compas_viewer import Viewer
import time
import cProfile
import pstats

from compas_fea2_vedo.viewer import ModelViewer


# ==============================================================================
# Set Up Profiling
# ==============================================================================
def profile_function(func, *args, **kwargs):
    """Profile a function and print a summary."""
    profiler = cProfile.Profile()
    profiler.enable()
    result = func(*args, **kwargs)
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.strip_dirs().sort_stats("cumtime").print_stats(
        10
    )  # Show top 10 slowest calls
    return result


# ==============================================================================
# Visualize in COMPAS Viewer
# ==============================================================================
def view_meshes(geo):
    viewer = Viewer()
    viewer.renderer.camera.scale = 100
    viewer.renderer.camera.position = [0, 0, 2000]
    viewer.renderer.camera.target = prt.bounding_box.to_mesh().centroid()
    viewer.renderer.camera.near = 10
    viewer.renderer.camera.far = 10000
    # viewer.scene.add(prt.discretized_boundary_mesh, show_faces=True, show_edges=True)
    viewer.scene.add(geo, show_faces=True, show_edges=True, opacity=1, show_points=True)
    viewer.show()


def merge_all_faces(mesh):
    """Merge all faces of a mesh into a single face.

    Parameters
    ----------
    mesh : compas.datastructures.Mesh
        The input mesh.

    Returns
    -------
    Mesh
        A new mesh with a single merged face.
    """
    boundary_vertices = mesh.vertices_on_boundaries()[0]
    bv_coords = mesh.vertices_attributes("xyz", keys=boundary_vertices)
    merged_mesh = Mesh.from_vertices_and_faces(
        bv_coords, [[i for i in range(len(bv_coords))]]
    )

    return merged_mesh


def merge_all_faces_delaunay(mesh):
    """Merge all faces of a mesh into a single Delaunay-triangulated mesh.

    Parameters
    ----------
    mesh : compas.datastructures.Mesh
        The input mesh.

    Returns
    -------
    Mesh
        A new mesh with Delaunay-triangulated faces.
    """
    from compas.geometry import delaunay_triangulation

    # Step 1: Get boundary vertices
    boundary_vertices = mesh.vertices_on_boundaries()[0]  # Extract first boundary loop
    bv_coords = mesh.vertices_attributes("xyz", keys=boundary_vertices)

    # Step 2: Compute Delaunay triangulation
    faces = delaunay_triangulation(
        [v[:2] for v in bv_coords]
    )  # Use only XY coordinates for 2D triangulation

    # Step 3: Create a new mesh with triangulated faces
    delaunay_mesh = Mesh.from_vertices_and_faces(bv_coords, faces)
    delaunay_mesh.unify_cycles()  # Ensure consistent face orientation

    return delaunay_mesh


units = units(system="SI_mm")

# ==============================================================================
# Define the data files
# ==============================================================================
HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "..", "00_data")
TEMP = os.path.join(HERE, "..", "..", "temp")

start_time = time.perf_counter()

# Set the backend implementation
# compas_fea2.set_backend("compas_fea2_opensees")
compas_fea2.set_backend("compas_fea2_castem")
print(f"Initialization took {time.perf_counter() - start_time:.4f} seconds")

mdl = Model(name="boxes")
mat_steel = Steel.S355(units=units)
sec = SolidSection(material=mat_steel)

prt = Part.from_step_file(
    step_file=os.path.join(DATA, "solids", "box.stp"),
    meshsize_max=200,
    section=sec,
)

mdl.add_part(prt)
print(f"Meshing took {time.perf_counter() - start_time:.4f} seconds")


planes = profile_function(
    prt.extract_clustered_planes, tol=1, angle_tol=2, verbose=False
)
submeshes = profile_function(
    prt.extract_submeshes,
    planes,
    tol=10,
    normal_tol=0.001,
    split=True,
)

submeshes = [merge_all_faces(submesh) for submesh in submeshes]


view_meshes(submeshes)
