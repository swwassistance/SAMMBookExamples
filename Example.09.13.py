###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
import scipy as sp

###############################################################################
# %% Example 9.13: Code for the implementation of the model of the 
#                  adsorption column 
###############################################################################
# Time
STARTTIME = 0                   # [d] Beginning
STOPTIME = [50, 100, 150, 200]  # [d] End of simulation       
DT = 0.1                        # [d] Time step

S = {}
q = {}
Sequ = {}
adsorbed = {}
for i in range(0, len(STOPTIME)):
    # Time span
    time = np.arange(STARTTIME, STOPTIME[i], DT)
    
    # Parameters: Process
    Q = 2000            # [m3 d-1] Influent
    S_in = 2            # [gDOC m-3] Adsorbable influent concentration
    kla = 5000          # [d-1] Rate constant for adsorption
    A = 10              # [m2] Cross section of the column
    eps = 0.4           # [-] Porosity of the carbon bed
    rho = 420000        # [gAC m-3reactor] Space density of the carbon
    qmax = 0.3          # [gDOC gAC-1] Maximum loading
    Ks = 10             # [gDOC m-3] Saturation concentration
    H = 3               # [m] Height of the column
    n = 30              # [-] Number of discrete nodes
    D = n*Q/(H*A*eps)   # [d-1] Transport rate per element
    
    # Parameters: Initial condition
    initS = np.zeros(n)     # [gDOC m-3] Concentration profile   
    initq = np.zeros(n)     # [gDOC gAC-1] Loading of the carbon
    initSequ = np.zeros(n)  # [gDOC m-3] Equilibrium concentration
    
    # Define ODE
    def model(var, t, param):
        S, q, Sequ = var
        Ks, qmax, kla, rho, D, S_in, eps = param
        
        Sequ = list(range(n))
        dqdt = list(range(n))
        dSdt = list(range(n))
        
        for i in range(n):
            # [gDOC m-3] Equilibrium concentration, Eq. 9.41
            Sequ[i] = q[i]*Ks/(qmax-q[i])                 
            dqdt[i] = kla*(S[i]-Sequ[i])/rho
            # Eq. 9.39 considering boundary condition
        dSdt[0] = D*(S_in-S[0])-kla*(S[0]-Sequ[0])/eps    
        dSdt[1:n] = D*(S[0:n-1]-S[1:n])-kla*(S[1:n]-Sequ[1:n])/eps   
        return dSdt, dqdt, Sequ
    
    # Solve ODE
    S[f'result{i+1}'], q[f'result{i+1}'], Sequ[f'result{i+1}'] = sh.sol_ode(
        model, var0=[initS, initq, initSequ], t=time, 
        param=[Ks, qmax, kla, rho, D, S_in, eps])          
 
###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 6.19, Fig. 9.14
fig = plt.figure('Example 9.13, Fig 9.14a')
plt.title('Example 9.13, Fig 9.14a')
for i in range(0, len(STOPTIME)):
    plt.plot([x/n*H for x in list(range(0, n))], S[f'result{i+1}'][-1, :],
             label=f't = {STOPTIME[i]} d')
plt.xlabel('Depth of adsorption column, z [m]')
plt.ylabel('Pollutant concentration S [gDOC m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

fig = plt.figure('Example 9.13, Fig 9.14b')
plt.title('Example 9.13, Fig 9.14b')
for i in range(0, len(STOPTIME)):
    plt.plot([x/n*H for x in list(range(0, n))], q[f'result{i+1}'][-1, :],
             label=f't = {STOPTIME[i]} d')
plt.xlabel('Depth of adsorption column, z [m]')
plt.ylabel('Load on activated carbon q [gDOC gAC$^{-1}$]')
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 6.19, Fig. 9.15

# Python version <= 3.10
i = len(STOPTIME)
fig = plt.figure('Example 9.13, Fig 9.15')
plt.title('Example 9.13, Fig 9.15')
S_in_vec = np.zeros(len(S[f'result{i}'][:, -1]))
S_in_vec[:] = 2
plt.plot(time, sp.integrate.cumulative_trapezoid(
    Q*(S_in_vec-S[f'result{i}'][:, -1])/1000, time, initial = 0), 
    label='$M_{DOC}$', color='black')
plt.ylabel('Sum of adsorbed pollutants $M_{DOC}$ [kgDOC]')
plt.xlabel('Time t [d]')
plt.ylim(0, 1000)
plt.yticks((0, 250, 500, 750, 1000))
plt.legend(loc=0, bbox_to_anchor=(1, 0.3))
plt.grid()
ax2 = plt.twinx()
plt.plot(time, S[f'result{i}'][:, -1], axes=ax2, label='S$_{out}$', 
         color='black', linestyle='--')
plt.ylabel('Effluent concentration S$_{out}$ [gDOC m$^{-3}$]')
plt.xlabel('Time t [d]')
plt.ylim(0, 2)
plt.yticks((0, 0.5, 1, 1.5, 2))
plt.legend(loc=1, bbox_to_anchor=(1, 0.2))
plt.show()

# # Python version > 3.10
# i = len(STOPTIME)
# fig, ax1 = plt.subplots(num ='Example 9.13, Fig 9.15',)
# ax1.set_title('Example 9.13, Fig 9.15')
# S_in_vec = np.zeros(len(S[f'result{i}'][:, -1]))
# S_in_vec[:] = 2
# ax1.plot(time, sp.integrate.cumulative_trapezoid(
#     Q*(S_in_vec-S[f'result{i}'][:, -1])/1000, time, initial = 0), 
#     label='$M_{DOC}$', color='black')
# ax1.set_ylabel('Sum of adsorbed pollutants $M_{DOC}$ [kgDOC]')
# ax1.set_xlabel('Time t [d]')
# ax1.set_ylim(0, 1000)
# ax1.set_yticks((0, 250, 500, 750, 1000))
# ax1.legend(loc=0, bbox_to_anchor=(1, 0.3))
# ax1.grid()
# ax2 = ax1.twinx()
# ax2.plot(time, S[f'result{i}'][:, -1], label='S$_{out}$', 
#          color='black', linestyle='--')
# ax2.set_ylabel('Effluent concentration S$_{out}$ [gDOC m$^{-3}$]')
# ax2.set_xlabel('Time t [d]')
# ax2.set_ylim(0, 2)
# ax2.set_yticks((0, 0.5, 1, 1.5, 2))
# ax2.legend(loc=1, bbox_to_anchor=(1, 0.2))
# plt.show()



