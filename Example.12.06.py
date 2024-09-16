###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 12.6: Sensitivity functions
###############################################################################
# This example serves as a skeleton only

# Parameters: Time
STARTTIME = 
STOPTIME = 
DT =        

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process


# Parameters: Initial conditions


# Define ODE
def model(var, t, param):
    = var
    = param
    
    return dCAdt

# SOLVE ODE
= sh.sol_ode(model, var0=[], t=time, param=[])

# Relative-relative sensitivtiy
p_sens = sh.sensitivity(par                 # Parameter p
                        xlabel,             # 'time' or 'length'
                        model,              # Model
                        var0=[],            # Initial condition
                        t=time,             # Time
                        param=[],           # Parameters of the model
                        param_var0=[],      # Parameters of var0
                        x_ind=-1)/par       # If two variables (X,S), 
                                            # default returns last