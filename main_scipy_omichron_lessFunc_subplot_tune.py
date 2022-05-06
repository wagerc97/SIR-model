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

def helperPlotEpsilonB(epsilon):
    g = np.linspace(-10, 10, 50)
    plt.plot(g, (epsilon(g)))
    plt.show()

def SIR_modelA(y, t, beta, gamma, epsilon):
    """
    Computes the derivative of y at t.
    :param y: result
    :param t: time
    :param beta: infection rate
    :param gamma: recovery rate
    :return:
    """
    S, I, R, = y
    dS_dt = -beta*S*I + epsilon*R
    dI_dt = beta*S*I - gamma*I
    dR_dt = gamma*I - epsilon*R

    return [dS_dt, dI_dt, dR_dt]

def SIR_modelB(y, t, beta, gamma, epsilon):
    """
    Computes the derivative of y at t.
    :param y: result
    :param t: time
    :param beta: infection rate
    :param gamma: recovery rate
    :return:
    """
    S, I, R, = y
    dS_dt = -beta*S*I + epsilon(t)*R
    dI_dt = beta*S*I - gamma*I
    dR_dt = gamma*I - epsilon(t)*R

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

def subplotSIR(solA, solB):
    # Plot the 3 phase lines for S, I and R
    # Model A
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10, 8])
    fig.suptitle("SIR-model A")
    ax1.plot(t, solA[:, 0], label="S(t) Susceptible", color="orange")
    ax1.plot(t, solA[:, 1], label="I(t) Infected", color="red")
    ax1.plot(t, solA[:, 2], label="R(t) Recovered", color="green")
    ax1.set_ylabel("Ratio of Population")
    ax1.legend(); ax1.grid();
    # Model B
    ax2.plot(t, solB[:, 0], label="S(t) Susceptible", color="orange")
    ax2.plot(t, solB[:, 1], label="I(t) Infected", color="red")
    ax2.plot(t, solB[:, 2], label="R(t) Recovered", color="green")
    ax2.set_ylabel("Ratio of Population")
    ax2.set_xlabel("Time in days")
    ax2.legend(); ax2.grid();
    plt.show()


if __name__ == '__main__':
    # Start the simulation with initial values
    print("=== Start calculation ===")

    # Initial conditions
    N = 1e6         # population size
    S_0 = 0.9       # fraction of Susceptible at t(0)
    I_0 = 0.1       # fraction of Infected
    R_0 = 0.0       # fraction of Recovered
    beta = 0.35     # infection rate
    gamma = 0.1     # recovery rate

    # Model A considers a constant decay in immunity over time
    epsilonA = 0.001

    # Model B considers the decay in immunity as a function over time
    # epsilonB_sigmoid is a sigmoid function with a damper
    # close to a jump function
    epsilonB_sigmoid = lambda t: (1 / (1 + np.exp(-t))) * 0.01
    # epsilonB_sigmoid is a gaussian function with a damper
    epsilonB_gauss = lambda t: None
    # epsilonB_sigmoid is a poisson function with a damper
    epsilonB_poisson = lambda t: None

    epsilonB = epsilonB_sigmoid # choice for epsilonB

    # Time
    t_0 = 0  # starting time
    dt = 1  # time step -> 1 day
    ndays = 150  # duration of simulation -> 5 months = 150 days
    t = np.arange(start=t_0, stop=ndays, step=dt)

    print(f"N = {int(N)}")
    print(f"S_0 = {S_0}")
    print(f"I_0 = {I_0}")
    print(f"R_0 = {R_0}")
    print(f"beta = {beta}")
    print(f"gamma = {gamma}")
    print(f"epsilonA = {epsilonA}")
    print(f"epsilonB = {epsilonB}")
    print(f"t_0 = {t_0}")
    print(f"timestep dt = {dt} day(s)")
    print(f"time = {ndays} days / {ndays/30} months")

    # Calculate Model for given parameters
    #solutionA = calculate(t, [t_0, ndays], [S_0, I_0, R_0], beta, gamma, epsilonA, SIR_modelA)
    #solutionB = calculate(t, [t_0, ndays], [S_0, I_0, R_0], beta, gamma, epsilonB, SIR_modelB)

    # Solve coupled system of ODEs.
    solutionA = np.array(solve_ivp(
        fun=lambda t, y: SIR_modelA(y, t, beta, gamma, epsilonA),  # system of ODEs
        y0=[S_0, I_0, R_0],     # initial condition
        t_span=[t_0, ndays],    # time points
        t_eval=t,               # when to store computed solutions
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix

    solutionB = np.array(solve_ivp(
        fun=lambda t, y: SIR_modelB(y, t, beta, gamma, epsilonB),  # system of ODEs
        y0=[S_0, I_0, R_0],     # initial condition
        t_span=[t_0, ndays],    # time points
        t_eval=t,               # when to store computed solutions
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix


    # Stdout Results
    print("\n---------------------------------------------")
    #print("The zenith of the infected wave is at", round(max(solutionA[:, 1]), 4))
    print("Plot is ready!")

    # Plot results
    plotSIR(solutionB)
    #subplotSIR(solutionA, solutionB)






"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                    ADDITIONAL INFORMATION 

‘RK45’ (default): Explicit Runge-Kutta method of order 5(4) [1]. 
The error is controlled assuming accuracy of the fourth-order method, 
but steps are taken using the fifth-order accurate formula (local extrapolation is done). 
A quartic interpolation polynomial is used for the dense output [2]. 
Can be applied in the complex domain.
 - Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp @method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""