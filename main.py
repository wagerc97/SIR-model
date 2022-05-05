"""""""""""""""""""""""""""""""""""""""""""""""""""""
SirModel.py
Script for 250070 SE Seminar Applied PDE (2022S)
Interpreter: Python 3.7
Author: Clemens Wager, BSc
"""""""""""""""""""""""""""""""""""""""""""""""""""""
# A Python script to compute SIR model of infection
# dynamics of the Covid-19 pandemic in Austria

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#4th-Order Runge-Kutta Method
#def RK4solver(params):


# Plot the 3 phase lines for S, I and R
#def plotSIR(params):


if __name__ == '__main__':
    # Start the simulation with initial values
    print("=== Start calculation ===")

    # Initial values
    N = 1e6 # population size
    S_0 = 0.9 # fraction of Susceptible at t(0)
    I_0 = 0.1 # fraction of Infected
    R_0 = 0.0 # fraction of Recovered
    a = 2.0 # infection rate
    b = 0.5 # recovery rate
    # Model A considers a constant decay in immunity over time
    c_A = 0.1
    # Model B considers a sigmoid function for the decay in immunity over time
    c_B = lambda t: 1 / 1 + np.exp(-t)
    # Time
    t_0 = 0 # starting time
    dt = 0.001 # time step
    duration = 5 # duration of simulation
