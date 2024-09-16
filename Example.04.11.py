###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import tqdm as tq

###############################################################################
# %% Example 4.11: Random walk as a method for integration
###############################################################################
# Please note that this simulation takes a long time.
# To make it faster, reduce the number of particles to n = int(100) 

# Parameters: Time
STARTTIME = 0           # [d] Start of the simulation
STOPTIME = 1000         # [d] End of the simulation
DT = 0.1                # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
xmin = 1                # [cm] Start of the bottle neck (along x-axis) 
                        #      in overall 3 cm wide bottle 
xmax = 2                # [cm] End of the bottle neck (along x-axis) 
                        #      in overall 3 cm wide bottle 
dB = 3                  # [cm] Diameter of the bottle 
hB = 15                 # [cm] Height of the bottle 
hF = 20                 # [cm] Height of the bottle including neck
D = 1                   # [cm2 d-1] Diffusion coefficient
Dx = np.sqrt(2*DT*D)    # [cm] Standard deviation of the single step 
                        #      in the random walk
n = int(100)            # [# particles] Number of particles
# [# particles] Number of particles in the bottle without neck
nF = int(n*dB*hB/(dB*hB+(xmax-xmin)*(hF-hB)))   

# Parameters: Initial condition
x = np.zeros((len(time), n))
y = np.zeros((len(time), n))
in_bottle = np.zeros((len(time), n))

# Starting position of the particles

# Particles in the bottle without neck
x[0][0:nF] = np.random.random(nF)*dB     
# Particles in the bottle neck                   
x[0][nF:n] = xmin + np.random.random((n-nF))*(xmax-xmin)
# Height of the particles in the bottle    
y[0][0:nF] = np.random.random(nF)*hB    
# Height of the particles in the bottle neck                    
y[0][nF:n] = hB + np.random.random((n-nF))*(hF-hB)          

for i in tq.tqdm(range(1, len(time))):
    for j in range(0, n):
        # New x coordinate after random step
        Delx = x[i-1][j] + Dx*np.random.normal(0, 1) 
        # New y coordinate after random step        
        Dely = y[i-1][j] + Dx*np.random.normal(0, 1)         
        # The next steps test that the particles do not move outside of the 
        # bottle (boundary condition)
        if Dely < hB:
            if Delx > 0 and Delx < dB:
                x[i][j] = Delx
            else:
                x[i][j] = x[i-1][j]
        else:
            if Delx > xmin and Delx < xmax:
                x[i][j] = Delx
            else:
                x[i][j] = x[i-1][j]
        
        if y[i-1][j] > hF:
            y[i][j] = y[i-1][j]
        else:
            if Dely < 0:
                y[i][j] = y[i-1][j]
            else: 
                if Dely < hB:
                    if Delx > 0 and Delx < dB:
                        y[i][j] = Dely
                    else:
                        y[i][j] = y[i-1][j]
                else:
                    if Delx > xmin and Delx < xmax:
                        y[i][j] = Dely
                    else:
                        y[i][j] = y[i-1][j]
    # Marks all particles which are still in the bottle
    in_bottle = np.where(y < hF, 1, 0)                      

# Counts the particles which are still in the bottle
in_bottle_sum = np.sum(in_bottle/n, axis=1)                 

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 4.11
fig = plt.figure('Example 4.11')
plt.title('Example 4.11')   
plt.scatter(x, y, s=0.1, color='black', label='In')
plt.scatter(np.where(y > hF, x, np.nan), np.where(y > hF, y, np.nan), s=0.1, 
            color='red', label='Out')
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.xlim(0, 3)
plt.ylim(0, 22)
plt.xticks([0, 1, 2, 3])
plt.yticks([0, 15, 20])
plt.legend(loc='center', fontsize=8, markerscale=10)
plt.gca().set_aspect('equal')
plt.show()       

# %% Plot Example 4.11, Fig. 4.8
fig = plt.figure('Example 4.11, Fig. 4.8')
plt.title('Example 4.11, Fig. 4.8')   
plt.plot([0, 200], [0.95, 0.95], color='red', linestyle='--')
plt.plot(time, in_bottle_sum, color ='black')
plt.text(205, 0.95, '5% particles which diffuse rapidly out of bottle neck', 
         va='center')
plt.xlabel('Time [d]')
plt.ylabel('Fraction of particles remaining [-]')
plt.xlim(0, 1000)
plt.ylim(0, 1)
plt.grid()
plt.show()
        
