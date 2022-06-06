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
from scipy.integrate import solve_ivp


def SIR_model(y, t, beta, gamma, omichron):
    """
    Computes the derivative of y at t.
    :param y: result
    :param t: time
    :param beta: infection rate
    :param gamma: recovery rate
    :return:
    """
    S, I, R, = y    # list of ODEs

    #dS_dt = -beta*S*I
    #dI_dt = beta*S*I - gamma*I
    #dR_dt = gamma*I
    dS_dt = -beta*S*I + omichron*R
    dI_dt = beta*S*I - gamma*I
    dR_dt = gamma*I - omichron*R

    return [dS_dt, dI_dt, dR_dt]



def plotSIR(solution):
    # Plot the 3 phase lines for S, I and R
    plt.figure(figsize=[10, 8])
    plt.plot(t, solution[:, 0], label="S(t) Susceptible", color="orange")
    plt.plot(t, solution[:, 1], label="I(t) Infected", color="red")
    plt.plot(t, solution[:, 2], label="R(t) Recovered", color="green")
    plt.legend(); plt.grid(); plt.title("SIR-model")
    plt.xlabel("Time"); plt.ylabel("Ratio of Population")
    plt.show()


if __name__ == '__main__':
    # Start the simulation with initial values
    print("=== Start calculation ===")

    # Model A considers a constant decay in immunity over time
    c_A = 0.01
    #c_A = lambda t: 0.1*t
    # Model B considers a sigmoid function for the decay in immunity over time
    c_B = lambda t: 1 / 1 + np.exp(-t)

    # Initial conditions
    N = 1e6  # population size
    S_0 = 0.9  # fraction of Susceptible at t(0)
    I_0 = 0.1  # fraction of Infected
    R_0 = 0.0  # fraction of Recovered
    a = 2.0  # infection rate
    b = 0.5  # recovery rate
    beta = 0.35  # infection rate
    gamma = 0.1  # recovery rate
    omichron = c_A

    # Time
    t_0 = 0  # starting time
    h = dt = 1  # time step -> 1 day
    ndays = 300  # duration of simulation -> 5 months = 150 days
    t = np.arange(start=t_0, stop=ndays, step=dt)

    print(f"N = {int(N)}")
    print(f"S_0 = {S_0}")
    print(f"I_0 = {I_0}")
    print(f"R_0 = {R_0}")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c_A = {c_A}")
    print(f"c_B = {c_B}")
    print(f"t_0 = {t_0}")
    print(f"timestep dt = {dt} day(s)")
    print(f"time = {ndays} days / {ndays/30} months")


    # Solve coupled system of ODEs.
    solution = np.array(solve_ivp(
        fun=lambda t, y: SIR_model(y, t, beta, gamma, omichron),  # system of ODEs
        y0=[S_0, I_0, R_0],     # initial condition
        t_span=[t_0, ndays],    # time points
        t_eval=t,               # when to store computed solutions
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix

    # Stdout Results
    print("\n---------------------------------------------")
    print("The zenith of the infected wave is at", round(max(solution[:, 1]), 4))
    print("Plot is ready!")

    # Plot results
    plotSIR(solution)



"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                    ADDITIONAL INFORMATION 

‘RK45’ (default): Explicit Runge-Kutta method of order 5(4) [1]. 
The error is controlled assuming accuracy of the fourth-order method, 
but steps are taken using the fifth-order accurate formula (local extrapolation is done). 
A quartic interpolation polynomial is used for the dense output [2]. 
Can be applied in the complex domain.
 - Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp @method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""