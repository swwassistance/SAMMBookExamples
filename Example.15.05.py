###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import math
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 15.5: Response of a fish population to a toxic spill
###############################################################################
# Parameters: Time
STARTTIME = 0           # [min] Beginning of the simulation
STOPTIME = 500          # [min] End of simulation 
DT = 1                  # [min] Time step

# Time span
tau = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
nfish = 1000        # [-] Number of individual fish 
nruns = 1000        # [-] Number of runs 
C = np.zeros((len(tau), nruns))   
E = np.zeros((len(tau), nruns))     
Damage = np.zeros((nfish, nruns)) 
FractionDamaged = np.zeros(nruns) 

# Random choice of tolerable exposures, describes variability of fish
Et = np.zeros(nfish)
Et[0:nfish] = np.random.normal(750, 50, nfish) 

for i in range(0, nruns):
    # Parameters: initial condition
    initC = np.random.uniform(5, 15) 
    # Initial exposure
    initE = 0 
    # [min-1] Random choice of recession constant for pollutant concentration
    k = np.random.normal(0.02, 0.002)
    
    # Define ODE
    def model(var, t, param):
        C, E = var
        # Time course of pollutant concentration
        dCdt = -k*C
        # Exposure is equal for all fish
        dEdt = C
        return dCdt, dEdt

    # Solve ODE
    C[:,i], E[:,i] = sh.sol_ode(model, var0=[initC, initE], t=tau, param=[])
    
    # Check whether damage occurs
    for j in range(0, nfish):
        if Et[j] > E[-1,i]:
            Damage[j,i] = 0
        else: 
            Damage[j,i] = 1
            
# Fraction of damaged fish             
FractionDamaged = np.mean(Damage, axis=0)

# [%] Cumulative distribution of damaged fish
sorted_FractionDamaged = np.sort(FractionDamaged)*100
cumulative_FractionDamaged = np.arange(1, len(FractionDamaged)+1) / len(FractionDamaged) 

# Find the values at 20%
setarr = np.zeros(nruns)
setpoint = 20
setarr[:] = setpoint
diff = np.abs(sorted_FractionDamaged-setarr)
min_index = np.argmin(diff)
res = cumulative_FractionDamaged[min_index]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 15.5, Fig. 15.3
plt.figure('Example 15.5, Fig. 15.3')
plt.title('Example 15.5, Fig. 15.3')
plt.plot(sorted_FractionDamaged, cumulative_FractionDamaged, 
         label='cumulative probability distribution')
plt.plot([setpoint, setpoint], [cumulative_FractionDamaged[0], res], 
         color = 'black', linestyle = '--')
plt.plot([sorted_FractionDamaged[0], setpoint], [res, res], 
         color = 'black', linestyle = '--')
plt.xlabel('Fraction of fish that will be damaged [%]')
plt.ylabel('Probability, that a given fraction \n will not be exceeded [-]')
plt.xlim(-1,100)
plt.ylim(0.65,1)
plt.grid()
plt.legend()
plt.show()

