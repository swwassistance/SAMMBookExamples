###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# %% Example 14.10: Simulation of an AR(1) process
# %% Example 14.11: Computation of autocorrelation function
###############################################################################
# Parameters: Time
STARTTIME = 0           # Beginning of computation
STOPTIME = 1000         # End of computation 
DT = 1                  # Time step of data series

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters
a1 = 0.9        # Regression weight α1
sigh = 1        # Standard deviation of white noise ση
N = int((STOPTIME-STARTTIME)/DT) # Number of elements in time series
order = 100 # Highest order of autocorrelation to be computed, < N/4

# Initial condition
x = np.zeros(N)  
eta = np.zeros(N)
cov = np.zeros([N, order])
cov_eta = np.zeros([N, order])
rho_data = np.zeros(order)
rho_data_eta = np.zeros(order)
rho_th = np.zeros(order)
rho_th[0] = a1**0
x[0] = 0 # Initializing time series, could be chosen randomly
eta[0] = 0
np.random.seed(8) # Seed = 8 makes it reproducible

for i in range(1, N):
    eta[i] = np.random.normal(0, sigh, 1)[0]
    x[i] = a1*x[i-1]+eta[i]  # Time series xi 
    
for i in range(0, order):
    rho_th[i] = a1**i

# k > 0 covariance, k = 0 variance of elements
for k in range(0, order): 
    for i in range(k, N): 
        cov[i, k] = cov[i-1, k]+x[i]*x[i-k] 
        cov_eta[i, k] = cov_eta[i-1, k]+eta[i]*eta[i-k] 
        
for k in range(0, order):
    rho_data[k] = N/(N-k)*cov[-1, k]/cov[-1, 0]
    rho_data_eta[k] = N/(N-k)*cov_eta[-1, k]/cov_eta[-1, 0]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.10, Fig. 14.14
plt.figure('Example 14.10, Fig. 14.14')
plt.title('Example 14.10, Fig. 14.14')
plt.plot(time, x, 'o', color='black', label='x$_i$', linestyle='-')
plt.xlabel('Time i')
plt.ylabel('x$_i$')
plt.xlim(0, 500)
plt.ylim(-6, 6)
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.10, Fig. 14.15
plt.figure('Example 14.10, Fig. 14.15')
plt.title('Example 14.10, Fig. 14.15')
plt.plot(time, eta, 'o', color='black', label='η$_i$ = x$_i$ - α$_i$x$_{i-1}$',
         linestyle='-')
plt.xlabel('Time i')
plt.ylabel('η$_i$')
plt.xlim(0, 500)
plt.ylim(-6, 6)
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.10, Fig. 14.16
plt.figure('Example 14.10, Fig. 14.16')
plt.title('Example 14.10, Fig. 14.16')
slope, intercept = np.round(np.polyfit(x[0:len(x)-1], x[1:len(x)], 1), 2)
if intercept > 0:
    plt.scatter(x[0:len(x)-1], x[1:len(x)], color='black', 
                label=f'AR(1) Model \nn= {N}\nRegression line: \n'
                f'x$_i$ = {slope}x$_{{i-1}}$ + {intercept}')
else:
    plt.scatter(x[0:len(x)-1], x[1:len(x)], color='black', 
                label=f'AR(1) Model \nn= {N}\nRegression line: \n'
                'x$_i$ = {slope}x$_{{i-1}}$ - {np.abs(intercept)}')
plt.xlabel('x$_{i-1}$')
plt.ylabel('x$_i$')
plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.10, Fig. 14.17
plt.figure('Example 14.10, Fig. 14.17')
plt.title('Example 14.10, Fig. 14.17')
slope, intercept = np.round(np.polyfit(x[0:len(x)-1], x[1:len(x)], 1), 2)
plt.scatter(eta[0:len(eta)-1], eta[1:len(eta)], color='black')
plt.xlabel('η$_{i-1}$')
plt.ylabel('η$_{i}$')
plt.xlim(-6, 6)
plt.ylim(-6, 6)
plt.grid()
plt.show()

# %% Plot Example 14.10, Fig. 14.18
plt.figure('Example 14.10, Fig. 14.18')
plt.title('Example 14.10, Fig. 14.18')
plt.plot(range(0, order), rho_data, 'o', color='black', label='Modelled ρ(k,x)',
         linestyle='-')
plt.plot(range(0, order), rho_th, color='black', label='Theoretical ρ(k,x)',
         linestyle='--')
plt.xlabel('Order, lag k')
plt.ylabel('ρ(k,x)')
plt.xlim(0, 30)
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.10, Fig. 14.19
plt.figure('Example 14.10, Fig. 14.19')
plt.title('Example 14.10, Fig. 14.19')
plt.plot(range(0, order), rho_data_eta, 'o', color='black', 
         label='Modelled ρ(k,η)', linestyle='-')
plt.xlabel('Order, lag k')
plt.ylabel('ρ(k,η)')
plt.xlim(0, 30)
plt.grid()
plt.legend()
plt.show()