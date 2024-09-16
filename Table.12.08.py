###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Table 12.8: Execution of a Monte Carlo simulation for the computation of 
# the uncertainty in the model prediction after Eq. (12.39), without correlation
###############################################################################
# Parameters: Time
STARTTIME = 0           # [h] Beginning of the simulation
STOPTIME = 2            # [h] End of the simulation
DT = 0.05               # [h] Time step of the model evaluation

# Parameters: Process
n = 10          # [-] Number of observations
C0_mu = 100     # [g m-3] Mean of normally distributed C0
C0_sig = 10     # [g m-3] Standard deviation of normally distributed C0
k_mu = 2        # [h-1] Mean of normally distributed k
k_sig = 0.5     # [h-1] Standard deviation of normally distributed k
# Correlation of C0 and k (0.95 if correlated, 0 if non-correlated)
corr = 0.95
runs = 10
# Standard deviatin of the residuals, Eq. 12.41
Res_sig = (n/(n-2)*(1-corr**2)*k_sig**2)**0.5 

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Variation
# [gm-3] Stochastic choice of C0 
C0_MC = np.random.normal(C0_mu, C0_sig, runs)
# [h-1] Stochastic choice of k, with correlation
k_MC = k_mu + corr*k_sig/C0_sig*(C0_MC-C0_mu) + np.random.normal(0, Res_sig, runs)

# Parameters: Initial condition
def var0(param_var0):
    initC = np.zeros(1) 
    C0 = param_var0
    initC = C0
    return initC

# Define ODE
def model(var, t, param):
    C = var
    k = param  
    dCdt = -k*C # Model equation
    return dCdt

# Solve ODE
sol = sh.sol_ode(model, var0(C0_mu), t=time, param=[k_mu])

# Monte Carlo
results, mean, stddev = sh.MonteCarlo(model, var0, t=time, param=[k_MC], 
                                      param_var0=[C0_MC], x_ind=0)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.08
plt.figure('Table 12.08, Fig. 12.21 or Fig. 12.24')
plt.title('Table 12.08, Fig. 12.21 or Fig. 12.24')
plt.plot(time, mean, 'black', linewidth=2)
plt.plot(time, mean+stddev, 'b--', linewidth=1)
plt.plot(time, mean-stddev, 'b--', linewidth=1)
plt.legend(['mean', 'standard deviation'], loc='upper right', 
           bbox_to_anchor=(1, 0.9))
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.show()

# %% Plot Example 12.08
fig = plt.figure('Table 12.08, Fig. 12.20')
plt.title('Table 12.08, Fig. 12.20')
for i in range(0, runs):
    plt.plot(time, results[i], color='black', linestyle='-') 
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.ylim(0, 120)
plt.grid()
plt.show()