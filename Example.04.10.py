###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# %% Example 4.10: Two-dimensional random walk 
###############################################################################
# Parameters: Time
STARTTIME = 0           # [s] Beginning of the simulation
STOPTIME = 10           # [s] End of the simulation
DT = 0.1                # [s] Duration of a time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
D = 2*10**(-3)          # [mm2 s-1] Diffusion coefficient 
Dx = np.sqrt(2*D*DT)    # [mm] Standard deviation of the local step in x 
Dy = Dx                 # [mm] Standard deviation of the local step in y

# Parameters: Initial condition
x = list(range(len(time)))
y = list(range(len(time)))
z = list(range(len(time)))

# Starting position
x[0] = 0                                    # [mm] Local coordinate x
y[0] = 0                                    # [mm] Local coordinate y

for i in range(1, len(time)):
    x[i] = x[i-1] + Dx*np.random.normal(0, 1) # [mm] Stochastic local step x
    y[i] = y[i-1] + Dy*np.random.normal(0, 1) # [mm] Stochastic local step y

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 4.10
fig = plt.figure('Example 4.10')
plt.title('Example 4.10')
plt.plot(x, y, color='black', linewidth=0.5, zorder=1)
plt.scatter(x, y, c=z, cmap='rainbow', s=100, zorder=2)
for i in range(0, len(time)):
    plt.text(x[i], y[i], z[i], ha='center', va='center', fontsize=8)  
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.grid()

# %% Plot Example 4.10, Fig. 4.7
n = 1000
distance_1d = list(range(len(time)))
distance_mean = list(range(len(time)))
distance_std = list(range(len(time)))
distance_std_th = list(range(len(time)))
distance_2d = np.zeros((n, len(time)))

for a in range(0, n):
    x[0] = 0
    y[0] = 0
    for b in range(1, len(time)):
        distance_1d[b] = Dx*np.random.normal(0, 1)
        distance_1d_cumsum = np.cumsum(distance_1d)
    distance_2d[a, :] = distance_1d_cumsum
distance_mean = np.mean(distance_2d, axis = 0)
distance_std = np.std(distance_2d, axis = 0)
distance_std_th = Dx*np.sqrt(time/DT)

fig = plt.figure('Example 4.10, Fig. 4.7')
plt.title('Example 4.10, Fig. 4.7')
plt.plot(time, distance_mean, color='black', label='Simulated', zorder=0)
plt.plot(time, np.zeros(len(time)), color='red',linestyle='--', 
         label='Theoretical', zorder=1)
plt.plot(time, distance_mean-distance_std, color='black', zorder=2)
plt.plot(time, -distance_std_th, color='red', linestyle='--', zorder=3)
plt.plot(time, distance_mean+distance_std, color='black', zorder=4)
plt.plot(time, distance_std_th, color='red', linestyle='--', zorder=5)
plt.text(5, 0.02, 'Mean', ha='center', va='center')  
plt.text(5, -0.11, 'Mean - std', ha='center', va='center')  
plt.text(5, 0.11, 'Mean + std', ha='center', va='center')  
plt.text(0.1, -0.175, 'n = 1000 Particles', ha='left', va='center')  
plt.xlabel('Time [s]')
plt.ylabel('Distance [mm]')
plt.grid()
plt.legend()
plt.show()


        
        