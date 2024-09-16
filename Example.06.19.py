###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 6.19: Numerical modeling of closed turbulent PFR
###############################################################################
# Parameters: Time
STARTTIME = 0   # [d] Beginning of the forward integration
STOPTIME = 3    # [d] End of the integration     
DT = 0.0002     # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = 25          # [-]       Number of partial reactors 
Vtot = 500      # [m3]      Total volume of the reactor
V = Vtot/n      # [m3]      Volume of a partial reactor, all of equal size
Q = 1000        # [m3 d-1]  Influent
R = 5000        # [m3 d-1]  Recirculation between partial reactors
Cm = 100        # [g m-3]   Average value of the inlet concentration
A = 50          # [g m-3]   Amplitude of the inlet concentration
f = 1           # [d-1]     Frequency of the variation of the inlet 
                #           concentration
k = 10          # [d-1]     Rate constant for degradation

# Parameters: Initial condition
initC = np.zeros(n)
initC[0:n] = 1  # [g m-3] Initial values for the concentrations

# Define ODE 
def model(var, t, param):
    C = var
    k, Q, V, R, Cm, A, f, n = param
    C_in = Cm+A*np.sin(2*np.pi*t/f)     # [g m-3] Example of a time-dependent 
                                        #         inlet concentration
    
    dCdt = list(range(n))
    # Balance for first reactor (closed for turbulence)
    dCdt[0] = (Q*C_in+R*C[1]-(Q+R)*C[0])/V-k*C[0] 
    
    # Balance for reactors 2 to n-1 (open)
    dCdt[1:n-1] = ((Q+R)*(C[0:n-2]-C[1:n-1])+R*(C[2:n]-C[1:n-1]))/V-k*C[1:n-1]
    
    # Balance for the last reactor (closed)
    dCdt[n-1] = ((Q+R)*(C[n-2]-C[n-1]))/V-k*C[n-1]
    
    return dCdt

# Solve ODE
C = sh.sol_ode(model, var0=initC, t=time, param=[k, Q, V, R, Cm, A, f, n])

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 6.19, Fig. 6.19
fig = plt.figure('Example 6.19, Fig 6.19')
plt.title('Example 6.19, Fig 6.19')
plt.plot(list(range(1,n+1)), C[-1][:], 'o', color='black', label='C(n)',
         linestyle='-')
plt.xlabel('Reactor number [n]')
plt.ylabel('C [g m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 6.19
fig = plt.figure('Example 6.19')
plt.title('Example 6.19')
plt.plot(time,C[:,-1], color='black', label='C$_n$(t)')
plt.xlabel('Time t [d]')
plt.ylabel('C$_n$ [g m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

