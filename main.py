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


# 4th-Order Runge-Kutta Method
# https://www.youtube.com/watch?v=6tUPGWDJghI -> Example on computer #TODO do on Friday
def RK4solver(params):
    # Time Update
    # t_i+1 = t_i + h

    # Intermediate steps
    # K_1 = f(t_i, S_i)
    # K_2 = f(t_i + h/2, S_i + 1/2 * K_1 * h)
    # K_3 = f(t_i + h/2, S_i + 1/2 * K_2 * h)
    # K_4 = f(t_i + h, S_i + K_3 * h)

    # Merge all steps
    # S_i+1 = S_i + h/6 * (K_1 + 2*k_2 + 3*K_3 + 4*K_4)

    pass


# Initiate Calculation
def initCalculation():
    pass


# Plot the 3 phase lines for S, I and R
def plotSIR(params):
    pass


if __name__ == '__main__':
    # Start the simulation with initial values
    print("=== Start calculation ===")
    # Initial values
    N = 1e6  # population size
    S_0 = 0.9  # fraction of Susceptible at t(0)
    I_0 = 0.1  # fraction of Infected
    R_0 = 0.0  # fraction of Recovered
    a = 2.0  # infection rate
    b = 0.5  # recovery rate
    # Model A considers a constant decay in immunity over time
    c_A = 0.1
    # Model B considers a sigmoid function for the decay in immunity over time
    c_B = lambda t: 1 / 1 + np.exp(-t)
    # Time
    t_0 = 0  # starting time
    h = 0.5  # time step -> 1 day
    duration = 150  # duration of simulation -> 5 months = 150 days

    print(f"N = {N}")
    print(f"S_0 = {S_0}")
    print(f"I_0 = {I_0}")
    print(f"R_0 = {R_0}")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c_A = {c_A}")
    print(f"c_B = {c_B(t_0)}")
    print(f"t_0 = {t_0}")
    print(f"h = {h} day(s)")
    print(f"duration = {duration} days / {duration / 30} months")
