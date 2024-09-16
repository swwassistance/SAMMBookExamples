###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Table 7.3: Implementation of a stochastic model of a turbulent plug-flow 
# reactor based on a random walk
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Beginning of the simulation
STOPTIME = 1            # [d] End of the simulation
DT = 0.001              # [d] Time step, must be small, see Eq. (7.44)

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
L = 10                  # [m] Length of the reactor
A = 1                   # [m2] Cross section of the reactor
Q = 20                  # [m3d-1] Flow rate
u = Q/A                 # [md-1] Flow velocity, advection
D = 5                   # [m2d-1] Dispersion coefficient
sig = (2*D*DT)**0.5     # [m] Standard deviation of the random step, Eq. (4.9)
nP = 10000              # [-] Number of particles

# Parameters: Initial condition
# [m] Initial location of particles without boundary conditions  
xnew = np.zeros((len(time), nP))
# [m] Initial location of particles with boundary conditions
x = np.zeros((len(time), nP))  
# For particles in the effluent
Effluent = np.zeros((len(time), nP))  
InEffluent = np.zeros(len(time))  
f = np.zeros(len(time))  

for i in range(1, len(time)):
    for ii in range(0, nP):
        # [m] New location of particles without boundary conditions
        xnew[i, ii] = x[i-1, ii] + np.random.normal(u*DT, sig)
        # Examination of the boundary conditions
        if xnew[i, ii] < 0:
            x[i, ii] = x[i-1, ii]
        elif xnew[i, ii] > L:
            x[i, ii] = L
        else:
            if x[i-1, ii] == L:
                x[i, ii] = L
            else:
                x[i, ii] = xnew[i, ii]
        # Marking the particles in the effluent
        if x[i, ii] == L:
            Effluent[i, ii] = 1
    # Counting the particles in the effluent
    InEffluent[i] = np.sum(Effluent[i, :], axis = 0) / nP

# Numeric derivative
f = (InEffluent - sh.delay(time, InEffluent, DT)) / DT

###############################################################################
# %% Plots
###############################################################################
# %% Plot Table 7.3
fig = plt.figure('Table 7.3a, Fig. 7.21')
plt.title('Table 7.3a, Fig. 7.21')
plt.plot(time, InEffluent, color='black', 
         label = f'Stochastic model computed with n = {nP}')
plt.xlabel('Space time $\\tau$ [d]')
plt.ylabel('Fraction of particles in the effluent, \n cum. RTD F($\\tau$) [-]')
plt.grid()
plt.show() 

# %% Plot Table 7.3
fig = plt.figure('Table 7.3b')
plt.title('Table 7.3b')
plt.plot(time, f, color='black')
plt.xlabel('Space time $\\tau$ [d]')
plt.ylabel('Numeric derivative [-]')
plt.grid()
plt.show()
