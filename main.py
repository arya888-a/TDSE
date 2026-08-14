import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#display in notebooks
from IPython.display import HTML, display

#Simulation Parameters
N = 1024                 # Number of spatial grid points
L = 200.0                # Length of the simulation box (-L/2 to L/2)
dx = L / (N - 1)         # Spatial step size
x = np.linspace(-L/2, L/2, N)

dt = 0.0005              # Time step (must be very small for Euler)
T_max = 30.0             # Final simulation time
n_steps = int(T_max / dt)
n_plot = 100             # Save a frame for animation every n_plot steps

print(f"Grid points: {N}, Spatial step: {dx:.4f}")
print(f"Time steps: {n_steps}, Time step: {dt:.6f}")

#Potential (V) and Initial State (Psi)

# A) Potential: Square barrier
V0 = 10.0                # Potential barrier height
barrier_width = 2.0      # Width of the barrier
barrier_center = 15.0    # Center of the barrier
V = np.zeros(N)
V[(x > barrier_center - barrier_width/2) & (x < barrier_center + barrier_width/2)] = V0

#Initial State: Gaussian Wave Packet
x0 = -30.0               # Initial position of the wave packet
p0 = 4.0                 # Initial momentum (velocity)
sigma = 3.0              # Width of the wave packet

# Create the wave packet
Psi_t = np.exp(-(x - x0)**2 / (4 * sigma**2)) * np.exp(1j * p0 * x)
Psi_t = Psi_t.astype(complex) # Ensure complex type

# Normalize the initial wave function
norm_initial = np.sum(np.abs(Psi_t)**2) * dx
Psi_t = Psi_t / np.sqrt(norm_initial)
print(f"Initial Norm: {np.sum(np.abs(Psi_t)**2) * dx:.6f}")

#Function to Apply the Hamiltonian (H * Psi)
#This function calculates the action of H on the wave function

def apply_H(Psi_in, V_pot, dx_val):
    # Calculate the second derivative (Laplacian) using Finite Difference
    # We use np.roll to get neighbors (Psi[j+1] and Psi[j-1])
    Psi_jp1 = np.roll(Psi_in, -1) # Psi[j+1]
    Psi_jm1 = np.roll(Psi_in, 1)  # Psi[j-1]
    
    d2Psi_dx2 = (Psi_jp1 - 2 * Psi_in + Psi_jm1) / (dx_val**2)
    
    # Apply hard-wall (Dirichlet) boundary conditions
    # This prevents the "wrap-around" from np.roll
    d2Psi_dx2[0] = 0
    d2Psi_dx2[-1] = 0
    
    # Calculate H * Psi
    H_Psi = -0.5 * d2Psi_dx2 + V_pot * Psi_in
    return H_Psi

#Time Evolution Loop (Forward Euler)

frames = []           # List to store animation frames
norms = []            # List to store the norm over time
times = []            # List to store time points

for i in range(n_steps):
    # Forward Euler update step
    H_Psi_t = apply_H(Psi_t, V, dx)
    Psi_t = Psi_t - 1j * dt * H_Psi_t
    
    # Enforce hard-wall boundary conditions (box)
    Psi_t[0] = 0
    Psi_t[-1] = 0
    
    # Store data for plotting
    if i % n_plot == 0:
        prob_density = np.abs(Psi_t)**2
        frames.append(prob_density)
        
        current_norm = np.sum(prob_density) * dx
        norms.append(current_norm)
        times.append(i * dt)
        
        # Print progress status
        if i % (n_steps // 10) == 0: 
            print(f"Step: {i}/{n_steps}, Time: {i*dt:.2f}, Norm: {current_norm:.6f}")

print("Simulation finished.")

#Plot of |Psi|^2 integral (Norm conservation)
print("\n--- Plotting Norm vs. Time ---")
plt.figure(figsize=(10, 5))
plt.plot(times, norms, 'r-')
plt.title("Norm Conservation Check (Forward Euler Method)", fontsize=16)
plt.xlabel("Time (t)", fontsize=14)
plt.ylabel("Norm ($\int |\Psi|^2 dx$)", fontsize=14)
# Start y-axis slightly below 1 to better show the increase
plt.ylim(0.99, max(norms) * 1.01) 
plt.grid(True)
plt.show()

#Animation of Time Evolution
print("--- Creating Animation ---")
fig, ax = plt.subplots(figsize=(12, 6))

ax.set_xlim(L/(-2), L/2)
ax.set_ylim(0, np.max(frames[0]) * 1.5) # Plot ceiling
ax.set_xlabel("Position (x)", fontsize=14)
ax.set_ylabel("Probability Density ($|\Psi|^2$)", fontsize=14)
ax.set_title("Wave Packet Time Evolution (Forward Euler)", fontsize=16)

#Plot the potential (scaled for visibility)
V_plot_scale = np.max(frames[0]) / (V0 * 1.5) if V0 > 0 else 0
ax.plot(x, V * V_plot_scale, 'k--', label=f"Potential (V0={V0})")
ax.legend(loc='upper right')

#The main line that will be animated
line, = ax.plot(x, frames[0], 'b-', lw=2, label="$|\Psi(x,t)|^2$")
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
norm_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)

def update_animation(frame_index):
    line.set_ydata(frames[frame_index])
    time_text.set_text(f"Time = {times[frame_index]:.2f}")
    norm_text.set_text(f"Norm = {norms[frame_index]:.4f}")
    
    # If instability causes the plot to explode, update the y-axis limit
    current_max = np.max(frames[frame_index])
    if current_max > ax.get_ylim()[1] * 0.9:
        ax.set_ylim(0, current_max * 1.2)
        
    return line, time_text, norm_text

# Create the animation
# (Note: In some environments, only the final plt.show() is needed)
ani = FuncAnimation(fig, update_animation, frames=len(frames), blit=True, interval=20)

# To display in a Jupyter Notebook
# print("Displaying animation... (this may take a moment)")
# display(HTML(ani.to_html5_video()))

# To show in a standard Python script
plt.show()

# To save the animation (requires ffmpeg)
# try:
#     ani.save('tdse_forward_euler.mp4', writer='ffmpeg', fps=30)
#     print("Animation saved successfully to tdse_forward_euler.mp4")
# except Exception as e:
#     print(f"Could not save animation. Do you have ffmpeg installed? Error: {e}")
