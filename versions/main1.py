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
import scipy.integrate


# https://www.youtube.com/watch?v=6tUPGWDJghI -> Example on computer
def Rk4SolverForOneOde(f, t_0, h, val_0):
    """
    4th-Order Runge-Kutta Method
    :param f: ODE
    :param t_0: old time
    :param h: time step
    :param val_0: old function value
    :return:
        t_1: updated time
        val_1: updated function value
    """
    # Intermediate steps
    # K_1 = f(t_i, S_i)
    K_1 = f(t_0, val_0)
    # K_2 = f(t_i + h/2, S_i + 1/2 * K_1 * h)
    K_2 = f(t_0 + h/2, val_0 + 0.5 * K_1 * h)
    # K_3 = f(t_i + h/2, S_i + 1/2 * K_2 * h)
    K_3 = f(t_0 + h/2, val_0 + 0.5 * K_2 * h)
    # K_4 = f(t_i + h, S_i + K_3 * h)
    K_4 = f(t_0 + h, val_0 + K_3 * h)

    # Update function value
    # S_i+1 = S_i + h/6 * (K_1 + 2*k_2 + 3*K_3 + 4*K_4)
    val_1 = val_0 + h/6 * (K_1 + 2*K_2 + 2*K_3 + K_4)

    # Time Update
    # t_i+1 = t_i + h
    t_1 = t_0 + h

    return t_1, val_1

# https://www.youtube.com/watch?v=qTUhLrwpiDs
def Rk4SolverForThreeOdes(f, t_0, h, val_0):
    """
    4th-Order Runge-Kutta Method
    :param f: ODE
    :param t_0: old time
    :param h: time step
    :param val_0: old function value
    :return:
        t_1: updated time
        val_1: updated function value
    """
    # solve with scipy.odeint: https://www.youtube.com/watch?v=wEvZmBXgxO0

    # Intermediate steps
    # K_1 = f(t_i, S_i)
    K_1 = f(t_0, val_0)
    # K_2 = f(t_i + h/2, S_i + 1/2 * K_1 * h)
    K_2 = f(t_0 + h/2, val_0 + 0.5 * K_1 * h)
    # K_3 = f(t_i + h/2, S_i + 1/2 * K_2 * h)
    K_3 = f(t_0 + h/2, val_0 + 0.5 * K_2 * h)
    # K_4 = f(t_i + h, S_i + K_3 * h)
    K_4 = f(t_0 + h, val_0 + K_3 * h)

    # Update function value
    # S_i+1 = S_i + h/6 * (K_1 + 2*k_2 + 3*K_3 + 4*K_4)
    val_1 = val_0 + h/6 * (K_1 + 2*K_2 + 2*K_3 + K_4)

    # Time Update
    # t_i+1 = t_i + h
    t_1 = t_0 + h

    return t_1, val_1


def initalCondition(S_0, I_0, R_0, t_0):
    S_1 = S_0
    I_1 = I_0
    R_1 = R_0
    t = t_0
    return S_1, I_1, R_1, t


def plotSIR(params):
    # Plot the 3 phase lines for S, I and R
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
    #c_A = 0.1
    c_A = lambda t: 0.1*t
    # Model B considers a sigmoid function for the decay in immunity over time
    c_B = lambda t: 1 / 1 + np.exp(-t)
    # Time
    t_0 = 0  # starting time
    h = dt = 1  # time step -> 1 day
    ndays = 150  # duration of simulation -> 5 months = 150 days

    # Define system of coupled nonlinear ODEs


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
