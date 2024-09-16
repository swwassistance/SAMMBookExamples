###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
from scipy.optimize import minimize

###############################################################################
# %% Example 14.7: Fitting a polynomial to another function
###############################################################################
# Parameters: Time
STARTTIME = 0           # [h] Beginning of simulation
STOPTIME = 5            # [h] End of simulation 
DT = 0.002              # [h] Time step

# Time span
tau = np.arange(STARTTIME, STOPTIME, DT)

# Parameters
n = 25          # Number of discretization steps
Vtot = 2500     # [m3] Volume of the reactor
V = Vtot/n      # [m3] Volume of a discrete element 
Q = 2500        # [m3 h-1] Influent and effluent flow
R = 4500        # [m3 h-1] Internal recirculation, turbulence
C_in = 0       

# Data import
data = sh.data_import('Example.14.07.Data.txt', ['Time', 'Nitrate'])
data_polynomial = data[60:172]
data_polynomial.reset_index(drop=True, inplace=True)

# Initial values of polynomial parameters
a = 1
b = 1
c = 1
d = 1

# Optimization
def objective_function1(params, fixed_params):
    a, b, c, d = params
    data = fixed_params
    chi2 = np.zeros(len(data.Nitrate))
    f = np.zeros(len(data.Nitrate))
    Pf = np.zeros(len(data.Nitrate))
    for i in range(0, len(data.Nitrate)):
        f[i] = a+b*data.Time[i]+c*data.Time[i]**2+d*data.Time[i]**3
        Pf[i] = data.Nitrate[i]
        chi2[i] = chi2[i-1]+(f[i]-Pf[i])**2
    chi2_sum = np.sum(chi2)
    return chi2_sum
op = minimize(objective_function1, [a, b, c, d], args=(data_polynomial)).x

# Optimization
def objective_function2(params, fixed_params):
    a, b, c, d = params
    data = fixed_params
    chi2 = np.zeros(len(data.Nitrate))
    f = np.zeros(len(data.Nitrate))
    Pf = np.zeros(len(data.Nitrate))
    for i in range(0, len(data.Nitrate)):
        f[i] = a+b*data.Time[i]+c*data.Time[i]**2+d*data.Time[i]**3
        Pf[i] = data.Nitrate[i]
        chi2[i] = chi2[i-1]+(f[i]-Pf[i])**2
    chi2_sum = np.sum(chi2)
    return chi2_sum, f, Pf

# Run at optimum
chi2_opt, f_opt, Pf_opt = objective_function2([op[0], op[1], op[2], op[3]], 
                                              data)

# Sign Changes
chi2sum_pol, f_pol, Pf_pol = objective_function2([op[0], op[1], op[2], op[3]], 
                                               data_polynomial)    
SignChange = np.zeros(len(data_polynomial))
SignChange[0] = 0

for i in range(1, len(data_polynomial.Nitrate)):
    if (f_pol[i]-Pf_pol[i])*(f_pol[i-1]-Pf_pol[i-1]) < 0:
        SignChange[i] = SignChange[i-1]+1
    else: 
        SignChange[i] = SignChange[i-1]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.7, Fig. 14.8
plt.figure('Example 14.7, Fig. 14.8')
plt.title('Example 14.7, Fig. 14.8')
plt.plot(data.Time, data.Nitrate, 'o', color='black', label='Measurements')
plt.plot(data_polynomial.Time, data_polynomial.Nitrate, 'o', color='red',
         label='Range of polynome fitting')
plt.plot(data.Time, f_opt, color='black', linestyle='--', label='Polynome')
plt.xlabel('Time [d]')
plt.ylabel('Nitrate, S [gN m$^{-3}$]')
plt.xlim(0.25, 1.5)
plt.ylim(2, 7)
plt.xticks([0.25, 0.5, 0.75, 1, 1.25, 1.5])
a = round(op[0], 1)
b = round(op[1], 1)
c = round(op[2], 1)
d = round(op[3], 1)
plt.text(0.3, 4.5, f'S = {a}{b}t+{c}t$^2${d}t$^3$', ha='left')
plt.grid()
plt.legend()
plt.show()

