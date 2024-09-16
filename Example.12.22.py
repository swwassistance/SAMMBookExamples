###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 12.22: MC simulation of the case study
###############################################################################
# Parameters: Time
STARTTIME = 0    # [h] Beginning of the simulation
STOPTIME = 45    # [h] End of the simulation
DT = 0.001       # [h] ime step of the model evaluation

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
runs = 100

# Parameters: Variation
# Stochastic Ssat
Ssat_mu = 9.85
Ssat_sig = 0.29
Ssat_MC = np.random.normal(Ssat_mu, Ssat_sig, runs)
Ssat_MC = sh.limit(Ssat_MC, 0) 
# Stochastic rO2max
rO2max_mu = -0.438
rO2max_sig = 0.016
rO2max_MC = np.random.normal(rO2max_mu, rO2max_sig, runs)
# Stochastic kla
kla_mu = 0.152
kla_sig = 0.0093
kla_MC = np.random.normal(kla_mu, kla_sig, runs) 
kla_MC = sh.limit(kla_MC, 0) 
# Stochastic S0
S0_mu = 2.03
S0_sig = 0.078
S0_MC = np.random.normal(S0_mu, S0_sig, runs)
S0_MC = sh.limit(S0_MC, 0)
# Stochastic KO
KO_mu = 0.996
KO_sig = 0.144
KO_MC = np.random.normal(KO_mu, KO_sig, runs)
KO_MC = sh.limit(KO_MC, 0)

# Parameters: Initial condition
def var0(param_var0):
    S0 = param_var0
    initS = np.zeros(1)
    initS = S0
    return initS

def model(var, t, param):
    S = var
    kla, Ssat, rO2max, KO = param
    # Reaction kinetics
    rO2 = rO2max*S/(KO+S)
    # Balance equation
    if t <= 15:
        dSdt = kla*(Ssat-S)+rO2
    else: 
        dSdt = rO2
    return dSdt

# Solve ODE
S = sh.sol_ode(model, var0([S0_mu]), t=time,
               param=[kla_mu, Ssat_mu, rO2max_mu, KO_mu])

# Monte Carlo
results, mean, stddev = sh.MonteCarlo(model, var0, t=time,
                                      param=[kla_MC, Ssat_MC, rO2max_MC, KO_MC],
                                      param_var0=[S0_MC])

# The second part with correlated parameters is missing, as P is unknown

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.22
plt.figure('Example 12.22')
plt.title('Example 12.22')
plt.plot(time, mean, 'black', linewidth=2)
plt.plot(time, mean+stddev, 'b--', linewidth=1)
plt.plot(time, mean-stddev, 'b--', linewidth=1)
plt.legend(['mean', 'standard deviation'], loc='upper right', 
           bbox_to_anchor=(1, 0.9))
plt.xlabel('Time [h]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.show()