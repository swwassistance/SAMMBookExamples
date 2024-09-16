###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 12.14: Computation of the error propagation after Eq. 12.37
###############################################################################
# Parameters: Time
STARTTIME = 0    # [min] Beginning of simulation
STOPTIME = 46    # [min] End of simulation
DT = 0.1         # [min] Time step, according to input

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
S_sat = 10      # [gO2 m-3] Estimated value for oxygen saturation
rO2 = -0.45     # [gO2 m-3 min-1] Estimated value for oxygen consumption
kla = 0.15      # [min-1] Estimated value for aeration coefficient
S0 = 2          # [gO2 m-3] Estimated value for initial oxygen concentration
tair = 15       # [min] Time at which the aeration is switched off
KO2 = 1         # [gO2 min-1] 
k = 5           # [-] Number of parameters

# Parameters: Covariance matrix (Table 12.6)
cov = np.array([[0.084, -0.003, -0.0026, 0.0113, 0.022], 
                [-0.003, 0.00026, 0.00008, -0.00014, -0.0022], 
                [-0.0026, 0.00008, 0.00009, -0.00044, -0.00056], 
                [0.0113, -0.00014, -0.00044, 0.0061, 0.00079], 
                [0.022, -0.0022, -0.00056, 0.00079, 0.021]])

# Parameters: Mean and standard error
mean = np.array([S_sat, rO2, kla, S0, KO2])
std = np.array([0.29, 0.016, 0.009, 0.078, 0.144])
string = np.array(['S$_{sat}$', 'r$_{O2}$', 'k$_l$a', 'S$_0$', 'K$_{O2}$'])

# Parameters: Initial condition
def var0(param_var0):
    S0 = param_var0
    initS = S0  # [gO2 m-3] Balance for oxygen concentration
    return initS       

# Define ODE
def model(var, t, param):
    S = var
    kla, S_sat, rO2, tair, KO2 = param
    S = sh.limit(S, 0) # Negative oxygen concentrations are not possible
    
    if t <= tair: # Switch off the aeration after 15 min
        on = 1
    else: on = 0
    
    dSdt = on*kla*(S_sat-S)+rO2*S/(KO2+S)
                  
    return dSdt

# Absolute-relative sensitivity
sens_ar = np.zeros((len(time), k))
sens_ar[:, 0] = sh.sensitivity(S_sat, 'time' ,model, var0, t=time,
                              param=[kla, S_sat, rO2, tair, KO2], param_var0=[S0])   
sens_ar[:, 1] = sh.sensitivity(rO2, 'time', model, var0, t=time,
                              param=[kla, S_sat, rO2, tair, KO2], param_var0=[S0])   
sens_ar[:, 2] = sh.sensitivity(kla, 'time', model, var0, t=time
                              ,param=[kla, S_sat, rO2, tair, KO2], param_var0=[S0])   
sens_ar[:, 3] = sh.sensitivity(S0, 'time', model, var0, t=time,
                              param=[kla, S_sat, rO2, tair, KO2], param_var0=[S0])   
sens_ar[:, 4] = sh.sensitivity(KO2, 'time', model, var0, t=time,
                              param=[kla, S_sat, rO2, tair, KO2], param_var0=[S0])   
    
# Relative-relative sensitivity
sens_rr = np.zeros((len(time), k))
for i in range(0, k):
    sens_rr[:, i] = sens_ar[:, i]/mean[i]

# Standard error of S, without covariance = Gaussian error propagation
A = np.zeros([len(time), k])
for i in range(0, k):
    A[:, i] = std[i]**2*sens_rr[:, i]**2
sigma_wo = np.sqrt(np.sum(A, axis=1))

# Standard error of S, with covariance
B = np.zeros([len(time), k])
sum1 = np.zeros([len(time), k])
C = np.zeros([len(time), k])
for i in range(0, k):
    for j in range(0, k):
        B[:, j] = sens_rr[:, j]*cov[i, j]
        sum1[:, i] = np.sum(B, axis=1)
    C[:, i] = sens_rr[:, i]*sum1[:, i]
sigma_w = np.sqrt(np.sum(C, axis=1))

# Solve ODE
initS = S0
S = np.array(sh.sol_ode(model, var0=[initS], t=time,
                        param=[kla, S_sat, rO2, tair, KO2]))
S = S.reshape(len(time))

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.14, Fig 12.14
plt.title('Example 12.14, Fig. 12.14')
plt.text(0.75, 6.8, 'Model E')
plt.xlabel('Time [min]')
plt.ylabel('Absolute relative sensitivity of S [g$_{O2}$ m$^{-3}$]')
plt.legend(string)
plt.grid(visible='bool')
plt.show()

# %% Plot Example 12.14, Fig 12.17
fig = plt.figure('Example 12.14, Fig. 12.17')
plt.title('Example 12.14, Fig. 12.17')
plt.plot(time, np.ones(len(time))*0.1, color='black', linestyle='-.',
         label='Standard error of individual measurements')
plt.plot(time, sigma_wo, color='black', linestyle='-',
         label='Without covariance = Gaussian error propagation')
plt.plot(time, sigma_w, color='black', linestyle='--',
         label='With covariance') 
plt.xlabel('Time [min]')
plt.ylabel('Standard deviation of S [g$_{O2}$ m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 12.14, Fig 12.18  
fig = plt.figure('Example 12.14, Fig. 12.18')
plt.title('Example 12.14, Fig. 12.18')
plt.plot(time, S, color='black', linestyle='-', label='S')
plt.plot(time, S+sigma_wo, color='black', linestyle='--',
         label='S+σ$_S$')
plt.plot(time, S-sigma_wo, color='black', linestyle='--',
         label='S-σ$_S$') 
plt.xlabel('Time [min]')
plt.ylabel('Standard deviation of S [g$_{O2}$ m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()