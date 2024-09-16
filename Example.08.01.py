###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 8.1: Effect of time of mixing on reactor performance
###############################################################################
# Parameters: Time
STARTTIME = 0           # [h]
STOPTIME = 5            # [h]
DT = 0.001              # [h]

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
orders = [0.5, 1, 2]            # [-] Reaction orders
constants = [20.03, 5, 0.534]   # [h-1] Rate constants 
n = 21                          # [-] Total number of compartments
Vtot = 1                        # [m3] Total volume
Q = 1                           # [m3 h-1] Influent flow rate 
C_in = 100                      # [g m-3] Influent concentration
       
eff = np.zeros((len(orders), n))
for i in range(len(orders)):
    m = orders[i]
    k = constants[i]
    
    # Parameters: Initial condition
    C = list(range(n))
    initC = np.zeros(n)
    V = np.zeros(n)

    for ii in range(n):
        iCSTR = ii
        # [g m-3] Initial concentration in all compartments
        initC[0:n] = 100      
        V[0:n] = Vtot/(2*(n-1))
        V[iCSTR] = Vtot/2
    
        # Define ODE
        def model(var, t, param):
            C = var
            Q, Vtot, C_in, m, k, iCSTR = param
            
            dCdt = list(range(n))   
            dCdt[0] = Q*(C_in - C[0])/V[0] - k*C[0]**m
            dCdt[1:n] = Q*(C[0:n-1] - C[1:n])/V[1:n] - k*C[1:n]**m
            return dCdt
    
        # Solve ODE
        C[ii] = sh.sol_ode(model, var0=initC, t=time, 
                           param=[Q, Vtot, C_in, m, k, iCSTR])

    eff[i][0:n] = [C[a][-1][-1] for a in range(n)]
    
###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 8.1, Fig. 8.5
fig = plt.figure('Example 8.1, Fig. 8.5')
plt.title('Example 8.1, Fig. 8.5')
for i in range(len(orders)):
    plt.plot(list(range(1,n+1)), eff[i][:], zorder=1)
    plt.scatter(list(range(1,n+1)), eff[i][:], 
                label = f'm = {orders[i]}, k = {constants[i]} h$^{{-1}}$', 
                zorder=2)
plt.text(1.2, 1.55, 'Early mixing',ha='left',va='bottom')
plt.text(20.8, 1.55, 'Late mixing',ha='right',va='bottom')
plt.text(11, 4.15, 'where:\nm = reaction order \nk = rate constant',
         ha='left', va='bottom')
plt.xlabel('Location of large CSTR, i$_{CSTR} [-]$')
plt.ylabel('Effluent pollutant concentration [M T$^{-1}$]')
plt.xlim(1, n)
plt.ylim(1.5, 5)
plt.xticks((1, 5, 9, 13, 17, 21))
plt.grid()
plt.legend()
plt.show()