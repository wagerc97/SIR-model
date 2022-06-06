"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
main.py
Script for 250070 SE Seminar Applied PDE (2022S)
Interpreter: Python 3.7
Author: Clemens Wager, BSc
----------------------------------------------------
# This is a Python script to simulate and visualize two SIR models.
# To model the infection dynamics of the Covid-19 pandemic in Austria.
# To test the effect on dynamics of decaying vaccination protection.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def SIR_modelA(y, t, beta, gamma, epsilon):
    """
    Computes the derivative of y at t. Callable by scipy.integrate.sove_ivp.
    :param y: result
    :param t: time
    :param beta: infection rate
    :param gamma: recovery rate
    :return:
    """
    S, I, R, = y
    #dS_dt = -beta*S*I + epsilon*R # epsilon term removable
    dS_dt = -beta*S*I
    dI_dt = beta*S*I - gamma*I
    #dR_dt = gamma*I - epsilon*R # epsilon term removable
    dR_dt = gamma*I

    return [dS_dt, dI_dt, dR_dt]

# epsilon is a function of time
def SIR_modelB(y, t, beta, gamma, epsilon):
    """
    Computes the derivative of y at t. Callable by scipy.integrate.sove_ivp.
    :param y: result
    :param t: time
    :param beta: infection rate
    :param gamma: recovery rate
    :return:
    """
    S, I, R, = y
    #print("type of t", type(t)) #numpy.float64
    #print("t =", t)
    #print("epsilon(t) =", epsilon(t))
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




def subplotsSIR1(solA, solB):
    # Plot the 3 phase lines for S, I and R
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,6))
    plt.setp((ax1, ax2), xticks=np.arange(0, ndays, 30), xticklabels=['Sep', 'Okt', 'Nov', 'Dec', 'Jan'],
             yticks=np.linspace(0, 1, 11))
    fig.suptitle("Comparison of SIR models\nModel A (decay rate = 0) vs Model B (decay rate as f(t))")
    # Model A
    ax1.plot(t, solA[:, 0], label="S(t) Susceptible", color="orange")
    ax1.plot(t, solA[:, 1], label="I(t) Infected", color="red")
    ax1.plot(t, solA[:, 2], label="R(t) Recovered", color="green")
    ax1.set_ylim([0.0, 1.0])
    ax1.set_xlim([0.0, ndays])
    ax1.set_ylabel("Ratio of Population")
    ax1.legend(); ax1.grid();
    # Model B
    ax2.plot(t, solB[:, 0], label="S(t) Susceptible", color="orange")
    ax2.plot(t, solB[:, 1], label="I(t) Infected", color="red")
    ax2.plot(t, solB[:, 2], label="R(t) Recovered", color="green")
    ax2.set_ylim([0.0, 1.0])
    ax2.set_xlim([0.0, ndays])
    ax2.set_ylabel("Ratio of Population")
    ax2.set_xlabel("Time")
    ax2.legend(); ax2.grid();
    plt.show()


def subplotsSIR(solA, solB):
    # Plot the 3 phase lines for S, I and R
    fig = plt.figure(figsize=(8,6))
    fig.suptitle("Comparison of SIR models\nModel A (decay rate = 0) vs Model B (decay rate as f(t))")
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)  # shared x axes
    plt.setp((ax1,ax2), xticks=np.arange(0,ndays,30))  # set xticks

     # Model A
    ax1.plot(t, solA[:, 0], label="S(t) Susceptible", color="orange")
    ax1.plot(t, solA[:, 1], label="I(t) Infected", color="red")
    ax1.plot(t, solA[:, 2], label="R(t) Recovered", color="green")
    ax1.set_ylim([0.0, 1.0])
    ax1.set_xlim([0.0, ndays])
    ax1.set_ylabel("Ratio of Population")
    ax1.legend(loc="upper right"); ax1.grid();
    # Model B
    plt.setp(ax1.get_xticklabels(), visible=False)  # dont show xlabels on ax1
    ax2.set_xticklabels(['Sep', 'Okt', 'Nov', 'Dec', 'Jan'])

    ax2.plot(t, solB[:, 0], label="S(t) Susceptible", color="orange")
    ax2.plot(t, solB[:, 1], label="I(t) Infected", color="red")
    ax2.plot(t, solB[:, 2], label="R(t) Recovered", color="green")
    ax2.set_ylim([0.0, 1.0])
    ax2.set_xlim([0.0, ndays])
    ax2.set_ylabel("Ratio of Population")
    ax2.set_xlabel("Time")
    ax2.legend(loc="upper right"); ax2.grid();
    plt.subplots_adjust(wspace=0.1, hspace=0.09)
    plt.show()



def helperPlotEpsilonB(epsilon, xstart, xstop):
    g = np.linspace(xstart, xstop, 50)
    plt.plot(g, (epsilon(g)))
    plt.title("Plot epsilon function for Model B\n(rate of protection decay)")
    plt.xticks(ticks=np.arange(xstart, xstop, 30), labels=['Sep', 'Okt', 'Nov', 'Dec', 'Jan'])
    plt.yticks(np.linspace(0, 1, 11)); plt.ylim(0, 1); plt.xlim(xstart, xstop)
    plt.grid(); plt.show()

def helperMap(t):
    # https://math.stackexchange.com/questions/377169/going-from-a-value-inside-1-1-to-a-value-in-another-range
    a = t_0
    b = ndays
    c = -3
    d = 0
    y = (t-a) * ((d-c)/(b-a)) + c
    return y


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

    I_0 = 1353/N       # fraction of Infected (1353 reported cases)

    R_0 = 0.58 + 0.0       # fraction of Recovered

    # 58% full vaccine protection, recovered?? vaccine age??
    S_0 = 1 - I_0 - R_0     # fraction of Susceptible at t(0)

    # k =  # contact rate
    # q =  # probability of an infection
    # D = 7 # duration of infectious state in days
    # beta = k * q * D # infection rate
    beta = 1.07 #0.35     # infection rate

    gamma = 10/ndays     # recovery rate (10 days)

    ##### Model A considers a constant decay in immunity over time #####
    epsilonA = 0.0

    ##### Model B considers the decay in immunity as a function over time #####
    # damps the function of epsilon
    damper = 0.05
    ### Constant rate of decay
    epsilonB_constant = lambda t: 1 / ndays * t * damper

    ### sigmoid function with a damper (sigmoid is similar to a jump function)
    epsilonB_sigmoid = lambda t: (1 / (1 + np.exp(-t)) ) * damper

    ### Half Normal Distributed values for decay rates
    epsilonB_halfNorm = lambda t: 1. / (np.sqrt(1 ** np.pi)) * np.exp(-1 * np.power(helperMap(t), 2.)) * damper

    ######### CHOICE #########
    epsilonB = epsilonB_halfNorm
   # helperPlotEpsilonB(epsilonB, xstart=0, xstop=150)


    print(f"N = {int(N)}")
    print(f"S_0 = {round(S_0, acc)}")
    print(f"I_0 = {round(I_0, acc)}")
    print(f"R_0 = {round(R_0, acc)}")
    print(f"damper = {damper}")
    print(f"beta = {round(beta, acc)}")
    print(f"gamma = {round(gamma, acc)}")
    print(f"epsilonA = {round(epsilonA, acc)}")
    #print(f"epsilonB = {round(np.mean(epsilonB(t)), acc)}")
    print(f"t_0 = {t_0}")
    print(f"timespan = {ndays} days / {ndays/30} months")
    print(f"timestep dt = {dt} day")


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
        t_eval=t,               # times when to store computed solution
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix


    # Stdout results
    print("\n---------------------------------------------")
    print("Day 0: 1st September, 2021")
    print(f"Day {ndays}: 31st January, 2021")
    print("Infected maximum in Model A:", round(max(solutionA[:, 1]), acc))
    print("Infected maximum in Model B:", round(max(solutionB[:, 1]), acc))
    print("Plot is ready!")

    # Plot results #
    #plotSIR(solutionB)
    subplotsSIR(solutionA, solutionB)





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                    ADDITIONAL INFORMATION 

‘RK45’ (default): Explicit Runge-Kutta method of order 5(4) [1]. 
The error is controlled assuming accuracy of the fourth-order method, 
but steps are taken using the fifth-order accurate formula (local extrapolation is done). 
A quartic interpolation polynomial is used for the dense output [2]. 
Can be applied in the complex domain.
 - Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp @method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""