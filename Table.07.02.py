###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Table 7.2: Implementation of a stochastic model of a cascade of stirred 
# tank reactors
###############################################################################
# Parameters: Time
STARTTIME = 0           # [T] Beginning of the simulation
STOPTIME = 1            # [T] End of the simulation
DT = 0.001              # [T] Time step, must be small such that p<1%

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
nR = 5                  # [-] Number of reactors 
Vtot = 1                # [m3] Total volume 
V = Vtot/nR             # [m3] Volume of a single reactor compartment
Q = 2                   # [m3/T] Flow rate
p = DT*Q/V              # [-] Probability that a particle leaves the reactor 
                        # compartment in the time step DT; 
                        # must be smaller than 0.01
nP = 10000              # Number of particles

# Parameters: Initial condition
loc = np.zeros((len(time),nP))   
                                 
loc[0,:] = 1
        
effluent = np.zeros((len(time),nP)) 
InEffluent = np.zeros((len(time)))

# Location of particles at time t (reactor comp.)                         
for i in range(1, len(time)):
    for ii in range(0, nP):
        rand = np.random.uniform(0, 1)
        if rand < p:
            loc[i, ii] = loc[i-1, ii] + 1
        else:
            loc[i, ii] = loc[i-1, ii]
        # Marking of the particles in the effluent
        if loc[i, ii] > nR:
            effluent[i, ii] = 1
        else:
            effluent[i, ii] = 0
    # Counting the particles already in the effluent 
    InEffluent[i] = np.sum(effluent[i, :]) / nP

# Numeric derivative, requires nP > 100,000
f = (InEffluent - sh.delay(time, InEffluent, DT)) / DT

###############################################################################
# %% Plots
###############################################################################
# %% Plot Table 7.2
fig = plt.figure('Table 7.2a')
plt.title('Table 7.2a')
plt.plot(time, InEffluent, color='black')
plt.xlabel('Time [T]')
plt.ylabel('Fraction in effluent [-]')
plt.grid()
plt.show()

# %% Plot Table 7.2
fig = plt.figure('Table 7.2b')
plt.title('Table 7.2b')
plt.plot(time, f, color='black')
plt.xlabel('Time [T]')
plt.ylabel('Numeric derivative [-]')
plt.grid()
plt.show()