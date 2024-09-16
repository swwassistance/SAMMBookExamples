###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 13.14: Modeling of a delay
###############################################################################
# Parameters: Time
STARTTIME = 0    # Beginning of the simulation
STOPTIME = 10    # End of the simulation
DT = 0.02        # Time step of the model evaluation

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = 6               # Order of the delay >= 1
THtot = 0.5         # Average delay
Tt = 0.25           # Dead time
TH = (THtot-Tt)/n   # Delay for one order

# Input signal Vin (= AR(1) model) as a function of time 
Vin = np.zeros(len(time))
Vin[0] = 2.5                        # Initial condition for Vin
for i in range(1, len(time)):       # Temporal change of Vin
    Vin[i] = 0.8 * Vin[i-1] + np.random.uniform(0, 1)
    
# Parameters: Initial condition
initV = np.zeros(n)

# Model
def model(var, t, param):
    V = var
    TH = param
    V_in = np.interp(t, time, Vin)
    dVdt = list(range(n))
    # Delay of the signal in nth order
    dVdt[0] = (V_in - V[0])/TH
    dVdt[1:n] = ((V[0:n-1] - V[1:n]))/TH
    return dVdt

# Solve ODE
V = sh.sol_ode(model, var0=initV, t=time, param=TH)

# Addition of the dead time, retarded signal
Vout = sh.delay(time, V, Tt)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.22
plt.figure('Example 13.14')
plt.title('Example 13.14, Fig. 13.12')
plt.plot(time, Vin, 'black', linestyle='dashed',label='Input')
plt.plot(time, Vout, 'black',label='Delayed output')
plt.xlabel('Time')
plt.ylabel('Input and delayed output signal')
plt.xlim(xmin=7, xmax=10)
plt.ylim(ymin=1,ymax=4)
plt.legend(loc='upper right')
plt.text(8.5,1.1, f'Dead time T$_t$ = {Tt}, θ$_{{tot}}$ = {THtot-Tt}, order n = {n}', 
         horizontalalignment='center', verticalalignment='center')
plt.show()
