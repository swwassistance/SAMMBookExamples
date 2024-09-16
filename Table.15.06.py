###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
from scipy.optimize import root

###############################################################################
#  %% Table 15.6: Stochastic simulation of the ozonation reactor
###############################################################################
# Inner loop only: uncer_runs = 1 and var_runs = 1000
# Both loops: uncer_runs = 1000 and var_runs = 1000 --> needs too much time

# Parameters: Time
STARTTIME = 0           # [d] Beginning of the simulation
STOPTIME = 0.2          # [d] Time to reach steady state
DT = 0.00025            # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)
    
# Parameters: Process
Vtots = [500, 600, 700, 800]    # [m3] Total reactor volume
uncer_runs = 1                  # [-] Runs uncertainty (= outer loop) 
var_runs = 1000                 # [-] Runs variation (= inner loop)
n = 6                           # [-] Number of reactors
Ccin = 1                        # [-] Relative influent concentration of cysts

res = np.zeros((len(Vtots), uncer_runs))
F = np.zeros((uncer_runs*var_runs))

for v in range(len(Vtots)):

    # Parameters: Process
    Vtot = Vtots[v]                  # [m3] Total reactor volume
    
    for u in range(0, uncer_runs):
        
        # Uncertainty
        fuQ = np.random.normal(1, 0.05)     # Normally distributed
        fuV = np.random.uniform(0.8, 1)     # Randomly distributed 
        fuSO31 = np.random.normal(1, 0.05)  # Measuring error, normally distributed
        fuSO34 = np.random.normal(1, 0.05)  # SO4 is independent of SO1
        fukO3 = 1                           # Uncertainty of kO3 is neglected
        fukD = 1                            # Uncertainty of kD is neglected
           
        # Variation
        fvSO31 = np.random.uniform(0.9, 1.1, var_runs)  # Evenly distributed
        fvSO34 = np.random.uniform(0.9, 1.1, var_runs)  # Evenly distributed
        fvkO3 = np.random.uniform(0.7, 1.3, var_runs)   # Randomly distributed
        fvkD = np.random.uniform(0.8, 1.2, var_runs)    # Randomly distributed
        
        Q = fuQ * 10000          # [m3d-1] Influent
        V = fuV * Vtot / n       # [m3] Active volume of a single reactor
        kO3 = 52 * fukO3 * fvkO3 # [d-1] Reaction constant for ozone decay at 5°C
        kD = 230 * fukD  * fvkD  # [m3g−1d−1] Reaction constant for disinfection at 5°C
        SO31 = fvSO31            # Controlled, remain constant
        SO34 = fvSO34            # Controlled, remain constant            
        
        # Parameters: Initial condition
        def var0(param_var0):
            SO31, SO34 = param_var0
            initSO3 = np.zeros(n)
            initCc = np.zeros(n)
            initSO3[0] = SO31
            initSO3[3] = SO34
            initSO3[1:3] = 1
            initSO3[4:6] = 1
            initCc[:] = 0.001   
            return initSO3, initCc
        
        # Define ODE
        def model(var, t, param):
            SO3, Cc = var
            Q, V, kO3, Ccin, kD = param   
            dSO3dt = list(range(n))
            dSO3dt[0] = 0
            for i in range(1, 3):
                dSO3dt[i] = Q * (SO3[i-1] - SO3[i]) / V - kO3 * SO3[i]
            dSO3dt[3] = 0
            for i in range(4, 6):
                dSO3dt[i] = Q * (SO3[i-1] - SO3[i]) / V - kO3 * SO3[i]        
            dCcdt = list(range(n))
            dCcdt[0] = Q * (Ccin - Cc[0]) / V - kD * Cc[0] * SO3[0]
            for i in range(1, 6):
                dCcdt[i] = Q * (Cc[i-1] - Cc[i]) / V - kD * Cc[i] * SO3[i]
            return dSO3dt, dCcdt
        
        # Monte Carlo
        results, mean, stddev = sh.MonteCarlo(model, var0, t=time, 
                                              param=[Q, V, kO3, Ccin, kD], 
                                              param_var0=[SO31, SO34])
        
        q = np.zeros(len(results))
        for i in range(0, len(results)):
            q[i] = results[i][-1]

        res[v,u] = np.percentile(q,95)

sort = np.sort(res)
F = np.array(range(uncer_runs))/((uncer_runs))

###############################################################################
# %% Plots
###############################################################################
# %% Plot Table 15.6
if uncer_runs > 1: 
    plt.figure('Table 15.6, Fig. 15.5')
    plt.title('Table 15.6, Fig. 15.5')
    for v in range(len(Vtots)):
        plt.plot(sort[v,:], F, label = f'V = {Vtots[v]} m$^3$')
    plt.xlabel('Residual fraction of remaining cysts C/C$_0$, 95% value')
    plt.ylabel('Cumulative frequency of the 95% value [-]')
    plt.grid()
    plt.legend()
    plt.show()