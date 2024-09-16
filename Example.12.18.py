###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 12.18: MC simulation with two correlated parameters, Option 1
###############################################################################
# Parameters: Time
STARTTIME = 0   # [h] Beginning of the simulation
STOPTIME = 2.05 # [h] End of the simulation
DT = 0.05       # [h] Time step of the model evaluation

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = 10          # [-] Number of observations
C0_mu = 100     # [g m-3] Mean of normally distributed C0
C0_sig = 33     # [g m-3] Standard deviation of normally distributed C0
k_mu = 2        # [h-1] Mean of normally distributed k
k_sig = 0.5     # [h-1] Standard deviation of normally distributed k
# Correlation of C0 and k (0.95 if correlated, 0 if non-correlated)
corr = 0.95
runs = 1000
# Standard deviatin of the residuals, Eq. 12.41
Res_sig = (n/(n-2)*(1-corr**2)*k_sig**2)**0.5 

# # Parameters: Initial condition
# def var0(param_var0):
#     initC = np.zeros(1)
#     C0 = param_var0
#     initC = C0
#     return initC

# Define ODE
def model(var, t, param):
    C = var
    k = param  
    dCdt = -k*C # Model equation
    return dCdt
    
# Monte Carlo
C_MC = np.zeros([len(time), runs])
for i in range(0, runs):
    # [g m-3] Stochastic variation of C0 
    C0 = np.random.normal(C0_mu, C0_sig)                                 
    # C = sh.limit(C0, 0)
    # [h-1] Stochastic variation of k
    k = k_mu+corr*k_sig/C0_sig*(C0-C0_mu)+np.random.normal(0, Res_sig)    
    # k = sh.limit(k, 0)
    C_MC[:, i] = sh.sol_ode(model, var0=C0, t=time, param=[k]).flatten()     

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.18
fig = plt.figure('Example 12.18a')
plt.title('Example 12.18a')
for i in range(0, runs):
    plt.plot(time, C_MC[:, i], color='black', linestyle='-') 
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.ylim(0, 120)
plt.grid()
plt.show()

# %% Plot Example 12.18
fig = plt.figure('Example 12.18b')
plt.title('Example 12.18b')
plt.plot(time, np.mean(C_MC[:], axis=1), color='black', linestyle='-', label='C')
plt.plot(time, np.mean(C_MC[:], axis=1)+np.std(C_MC[:], axis=1), color='black', 
         linestyle='--', label='C$\pm$σ$_C$')
plt.plot(time, np.mean(C_MC[:], axis=1)-np.std(C_MC[:], axis=1), color='black', 
         linestyle='--')
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

###############################################################################
# %% Example 12.18: MC simulation with two correlated parameters, Option 2
###############################################################################
# Parameters: Time
STARTTIME = 0   # [h] Beginning of the simulation
STOPTIME = 2.05 # [h] End of the simulation
DT = 0.05       # [h] Time step of the model evaluation

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = 10          # [-] Number of observations
C0_mu = 100     # [g m-3] Mean of normally distributed C0
C0_sig = 33     # [g m-3] Standard deviation of normally distributed C0
k_mu = 2        # [h-1] Mean of normally distributed k
k_sig = 0.5     # [h-1] Standard deviation of normally distributed k
# Correlation of C0 and k (0.95 if correlated, 0 if non-correlated)
corr = 0.95     
runs = 1000
# Standard deviatin of the residuals, Eq. 12.41
Res_sig = (n/(n-2)*(1-corr**2)*k_sig**2)**0.5 

# Parameters: Variation
C0_MC = np.random.normal(C0_mu, C0_sig, runs) 
k_MC = k_mu + corr*k_sig/C0_sig*(C0_MC-C0_mu)+np.random.normal(0, Res_sig, runs)

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
C = sh.sol_ode(model, var0(C0_mu), t=time, param=[k_mu])

# Monte Carlo
results, mean, stddev = sh.MonteCarlo(model, var0, t=time, param=[k_MC], 
                                      param_var0=[C0_MC], x_ind=0)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.18
plt.figure('Example 12.18b')
plt.title('Example 12.18b')
plt.plot(time, mean, 'black', linewidth=2)
plt.plot(time, mean+stddev, 'b--', linewidth=1)
plt.plot(time, mean-stddev, 'b--', linewidth=1)
plt.legend(['mean', 'standard deviation'], loc='upper right', 
           bbox_to_anchor=(1, 0.9))
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
