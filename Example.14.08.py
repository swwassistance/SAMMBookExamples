###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
import scipy as sp
from scipy.interpolate import interp1d

###############################################################################
# %% Example 14.8: Decomposition of a function based on Fourier transformation
###############################################################################
# Parameters: Time
STARTTIME = 0           # Beginning of integration
STOPTIME = 4.02         # Twice the length of the time series
DT = 1/50               # Time step
T = 2                   # Length of the time series

# Time spans
time = np.arange(STARTTIME, STOPTIME, DT)

# Time series example, see Fig. 14.13
time_meas = np.array([0, 0.5, 0.7, 0.8, 1.3, 1.5, 1.6, 2])
f_meas = np.array([25, 35, 7.5, 27.5, 30, 45, 25, 25])

# Parameters
n = 25          # Number of coefficients
m = 9           # Degree of the Fourier polynomial

Vtot = 2500     # [m3] Volume of the reactor
V = Vtot/n      # [m3] Volume of a discrete element
Q = 2500        # [m3 h-1] Influent and effluent flow
R = 4500        # [m3 h-1] Internal recirculation, turbulence
C_in = 0

# Parameters: Initial condition
inita = np.zeros(n)
initb = np.zeros(n)

# Define ODE
def model(var, t, param):
    a, b = var
    T = param
    dadt = list(range(n))
    dbdt = list(range(n))
    for i in range(n):
        if t < T:
            f = interp1d(time_meas, f_meas, kind='linear', 
                         fill_value='extrapolate')
            dadt[i] = 2*f(t)*np.cos(2*np.pi*i*t/T)/T
            dbdt[i] = 2*f(t)*np.sin(2*np.pi*i*t/T)/T
        else:
            dadt[i] = 0
            dbdt[i] = 0
    return dadt, dbdt

# Solve ODE
a, b = sh.sol_ode(model, var0=[inita, initb], t=time, param=[T])

c = np.zeros([len(time), n])
c[:, 0] = a[:, 0]/2
for i in range(1, n):
    c[:, 1:n] = np.sqrt(a[:, 1:n]**2+b[:, 1:n]**2)
    
g = np.zeros([len(time), m])
for i in range(len(time)):
    for im in range(1, m):
        if time[i] < T:
            g[i, 0] = 0
            g[i, im] = 0
        else: 
            g[i, 0] = c[i, 0]
            g[i, im] = g[i, im-1] + a[i, im]*np.cos(im*time[i]*2*np.pi/T)
            +b[i, im]*np.sin(im*time[i]*2*np.pi/T)
         
FourierPolynomial = g[:, -1]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.8, Fig. 14.13
plt.figure('Example 14.8, Fig. 14.13')
plt.title('Example 14.8, Fig. 14.13')
plt.plot(time_meas, f_meas, label='Original Function', color='black', 
         linestyle='--')
plt.plot([2, 2.5, 2.7, 2.8, 3.3, 3.5, 3.6, 4], f_meas, color='black',
         linestyle='--')
plt.plot(time, FourierPolynomial, label=f'Fitted function \n'
         'Fourier series (n = {m})', color = 'black', linestyle='-')
plt.xlabel('Time [d]')
plt.ylabel('Concentration [g m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()
