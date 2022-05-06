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

# epsilon is a function of time
def SIR_modelB1(y, t, beta, gamma, epsilon):
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
# epsilon is 1/150, so protection wares of after 5 months
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
    dS_dt = -beta*S*I + epsilon*R
    dI_dt = beta*S*I - gamma*I
    dR_dt = gamma*I - epsilon*R

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

    # Time
    t_0 = 0  # starting time
    dt = 1  # time step -> 1 day
    ndays = 150  # duration of simulation -> 5 months = 150 days
    t = np.arange(start=t_0, stop=ndays, step=dt)

    # Initial conditions
    N = 9e6         # population size
    acc = 3         # accuracy, number of decimal places in results
    S_0 = 0.9       # fraction of Susceptible at t(0)
    I_0 = 1353/N    # fraction of Infected
    R_0 = 0.0       # fraction of Recovered
    beta = 0.35     # infection rate
        #k =  # contact rate
        #q =  # probability of an infection
        #D = 7 # duration of infectious state in days
        #beta = k * q * D # infection rate
    gamma = 10/ndays     # recovery rate (10 days)

    # Model A considers a constant decay in immunity over time
    epsilonA = 0.0

    # Model B considers the decay in immunity as a function over time
    # epsilonB_sigmoid is a sigmoid function with a damper
    # sigmoid is similar to a jump function
    epsilonB_sigmoid = lambda t: (1 / (1 + np.exp(-t))) * 0.02
    # Maybe try constant, Gaussian, Poisson, ... decay rate
    epsilonB = 1/150

    print(f"N = {int(N)}")
    print(f"S_0 = {S_0}")
    print(f"I_0 = {I_0}")
    print(f"R_0 = {R_0}")
    print(f"beta = {beta}")
    print(f"gamma = {gamma}")
    print(f"epsilonA = {epsilonA}")
    # print(f"epsilonB = {round(np.mean(epsilonB(t)), acc)}")
    print(f"epsilonB = {round(epsilonB, acc)}")
    print(f"t_0 = {t_0}")
    print(f"timestep dt = {dt} day(s)")
    print(f"time = {ndays} days / {ndays/30} months")

    # Solve coupled system of ODEs using RK4
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

    # Scale results on Austrian population
    scaled = False
    if scaled:
        solutionA = solutionA*N
        solutionB = solutionB*N
        acc = 0

    # Stdout results
    print("\n---------------------------------------------")
    print("Results scaled:", scaled)
    print("Zenith of  infected wave in Model A is at", round(max(solutionA[:, 1]), acc))
    print("Zenith of  infected wave in Model B is at", round(max(solutionB[:, 1]), acc))
    print("Plot is ready!")

    # Plot results
    #plotSIR(solutionB)
    subplotSIR(solutionA, solutionB)






"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                    ADDITIONAL INFORMATION 

‘RK45’ (default): Explicit Runge-Kutta method of order 5(4) [1]. 
The error is controlled assuming accuracy of the fourth-order method, 
but steps are taken using the fifth-order accurate formula (local extrapolation is done). 
A quartic interpolation polynomial is used for the dense output [2]. 
Can be applied in the complex domain.
 - Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp @method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""