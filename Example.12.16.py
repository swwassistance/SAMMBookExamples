###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 12.16: Accident with toxic materials
###############################################################################
# Parameters: Time
STARTTIME = 0       # [min] Beginning of the simulation
STOPTIME =  1000    # [min] End of simulation
DT = 1              # [min] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
runs = 1000
Et_m = 750      # [g min m-3] Mean of initial pollutant concentration
Et_sig = 50     # [g min m-3] Standard deviation of initial pollutant 
                #             concentration
k_m = 0.02      # [min-1]     Mean of recession constant for pollutant 
                #             concentration
k_sig = 0.002   # [min-1]     Standard deviation of recession constant for 
                #             pollutant concentration

# Parameters: Initial condition
damage = list(range(runs))  # Check whether damage occurs
initE = 0                   # Exposure of the individual

# Model
for i in range(runs):
    # [g min m-3] Random choice of tolerable exposure
    Et = np.random.normal(Et_m, Et_sig)   
    # [g m-3] Random choice of initial pollutant concentration with limits 
    # 5 and 15 g m-3
    C0 = np.random.uniform(5, 15)            
    # [min-1] Random choice of recession constant for pollution concentration 
    k = np.random.normal(k_m, k_sig)         

    # Define ODE
    def model(var, t, param):
        C = var
        C0, k = param  
        C = C0*np.exp(-k*t)
        dEdt = C
        return dEdt

    # Solve ODE
    E = sh.sol_ode(model, var0=[initE], t=time, param=[C0, k])

    if E[-1] > Et:
        damage[i] = 1
    else: 
        damage[i] = 0
    
damages = np.sum(damage)   
damage_prob = damages/runs 
print(round(damage_prob*100, 2),'%')
