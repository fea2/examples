from compas_fea2.model import Steel, AngleSection, RectangularSection, ISection

# # Define the section
# section = AngleSection(w=100, h=200, t1=10, t2=10, material=Steel.S355())
# section = RectangularSection(w=100, h=200, material=Steel.S355())
section = ISection.HEA280(material=Steel.S355())
section.shape.plot()

# Define the loads
N = -100e3
Mx = 10e6
My = 10e4
Vx = 15e3
Vy = 10e3

# # Compute the stress at a point
# sigma, tau_x, tau_y = isection.compute_stress(N, Mx, My, Vx, Vy, 0, 0)
sigma, tau_x, tau_y = section.compute_stress(N, Mx, My, Vx, Vy, 0, 0)
print("Normal stress: ", sigma)
print("Shear stress in x: ", tau_x)
print("Shear stress in y: ", tau_y)

section.plot_stress_distribution(
    N=N, Mx=Mx, My=My, Vx=Vx, Vy=Vy, nx=100, ny=100, show_tau=False
)
