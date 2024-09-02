# SAMM Book Examples
The [Chairs of Urban Water Management](https://sww.ifu.ethz.ch/) translated the examples in the [course book](https://link.springer.com/book/10.1007/978-3-540-77278-1) of the lecture "Systems Analysis and Mathematical Modeling in Urban Water Management (SAMM)" from Berkeley Madonna to Python.

## Table of Contents
[Requirements](#requirements)  
[Overview](#overview)  
[Contact](#contact)  

## Requirements
We recommend to download: 
- [Spyder](https://www.spyder-ide.org/)
- [Python 3.10.7](https://www.python.org/downloads/)
- [pip package manager](https://bootstrap.pypa.io/get-pip.py): Save the document as get-pip.py in a folder of your choice. Open the folder, type cmd into the location bar at the top of the window and tap enter. The command prompt window is opened in which you run ```python get-pip.py```. pip is now installed on your system.

After this, open a command window and install the required packages by entering ```pip install matplotlib numpy pandas scipy sammhelper odeintw tqdm```.

## Overview

### Chapter 4: Transport Processes
- Example 4.10: Two-dimensional random walk
- Example 4.11: Random walk as a method for integration

### Chapter 6: Ideal Reactors
- Example 6.17: Implementation of a turbulent PFR
- Example 6.19: Numerical modeling of closed turbulent PFR
- Example 6.23: Implementation of an SBR

### Chapter 7: Hydraulic Residence Time Distribution
- Example 7.6: Convolution
- Example 7.11: Simulation of the RTD of a cascade of CSTRs for Figure 7.10
- Example 7.15: Simulation of the RTD of a closed turbulent PFR
- Table 7.2: Implementation of a stochastic model of a cascade of stirred tank reactors
- Table 7.3: Implementation of a stochastic model of a turbulent plug-flow reactor based on a random walk

### Chapter 8: Modeling of Real Reactors
- Example 8.1: Effect of time of mixing on reactor performance
- Example 8.3 with Data 8.3: Expected value and variance of the residence time distribution from data
- Example 8.5 with Data 8.5: Code for the direct identification of the model parameters

### Chapter 9: Heterogeneous Systems
- Example 9.6: Computation of the concentration profiles in an activated sludge floc 
- Example 9.13: Code for the implementation of the model of the adsorption column

### Chapter 11: Measurement and Measurement Uncertainty
- Example 11.2: Log-normal distribution of a random variable

### Chapter 12: Parameter Identification, Sensitivity and Error Propagation
- Example 12.4: Implementation of Eq. 12.11
- Example 12.6: Sensitivity functions
- Example 12.14: Computation of the error propagation after Eq. 12.37
- Example 12.15: Computation of π with stochastic simulation
- Table 12.8: Execution of a Monte Carlo simulation for the computation of the uncertainty in the model prediction after Eq. (12.39), without correlation
- Example 12.16: Accident with toxic materials
- Example 12.18: MC simulation with two correlated parameters
- Example 12.22: MC simulation of the case study

### Chapter 13: Process Control Engineering
- Example 13.14: Modeling of a delay
- Example 13.21: Implementation of a two-position controller (n = 1)

### Chapter 14: Time Series Analysis
- Example 14.2 with Data 14.2: Test for stationarity of a time series
- Example 14.3 with Data 14.3: Computation of the arithmetic moving average
- Example 14.5 with Data 14.5: Computation of the geometric moving average
- Example 14.7 with Data 14.7: Fitting a polynomial to another function
- Example 14.8: Decomposition of a function based on Fourier transformation
- Example 14.10 and 14.11: Simulation of an AR(1) process and computation of autocorrelation function
- Example 14.13 with Data 14.13: Analysis of the time series in Fig. 14.21
- Example 14.15 and 14.16 with Data 14.15: Code for the analysis of the time series

### Chapter 15: Design under Uncertainty
- Example 15.5: Response of a fish population to a toxic spill
- Example 15.13: Uncertain definition of distributions may lead to unforeseeable results
- Table 15.4: Deterministic simulation of the ozonation and disinfection reactor
- Table 15.6: Stochastic simulation of the ozonation reactor

## Contact
In case of questions or feedback, feel free to contact [sww@ifu.baug.ethz.ch](sww@eifu.baug.ethz.ch).
