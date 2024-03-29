"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
main.py
Script for 250070 SE Seminar Applied PDE (2022S)
Interpreter: Python 3.9
Author: Clemens Wager, BSc
Last revisited: July 10th, 2022
---------------------------------------------------------------------------------------------------------
# This is a Python script to simulate and visualize two SIR models.
# To model the infection dynamics of the Covid-19 pandemic in Austria.
# To test the effect of decaying vaccination protection on pandemic dynamics.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_ivp


def SIR_modelA(y, t, beta, gamma, epsilon, immunity_level):
    """
    Computes the derivative of y at t. Callable by scipy.integrate.solve_ivp.
    Model A looks at constant vaccination protection without decay over time.
    :param y: result
    :param t: time
    :param beta: transmission rate
    :param gamma: recovery rate
    :param epsilon: constant float
    :param immunity_level: ratio of population that gains immunity [0-1]
    :return:
    """
    S, I, R, = y
    dS_dt = -beta * S * I + epsilon * R + (1 - immunity_level) * gamma * I
    dI_dt = beta * S * I - gamma * I 
    dR_dt = immunity_level * gamma * I - epsilon * R

    return [dS_dt, dI_dt, dR_dt]


def SIR_modelB(y, t, beta, gamma, epsilon, immunity_level):
    """
    Computes the derivative of y at t. Callable by scipy.integrate.solve_ivp.
    Model B looks at time dependent vaccination protection.
    :param y: result
    :param t: time
    :param beta: transmission rate
    :param gamma: recovery rate
    :param epsilon: time dependent function
    :param immunity_level: ratio of population that gains immunity [0-1]
    :return:
    """
    S, I, R, = y
    dS_dt = -beta * S * I + epsilon(t) * R + (1 - immunity_level) * gamma * I
    dI_dt = beta * S * I - gamma * I
    dR_dt = immunity_level * gamma * I - epsilon(t) * R

    return [dS_dt, dI_dt, dR_dt]



def subplotsSIR(solA, solB):
    # Plot the 3 phase lines for S, I and R
    fig = plt.figure(figsize=(8,6))
    fig.suptitle("Comparison of SIR models\nModel A (decay rate = 0) vs Model B (decay rate as f(t))")
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)  # shared x axes
    plt.setp((ax1,ax2), xticks=np.arange(0,ndays,30))  # set xticks

     # Model A
    ax1.plot(t, solA[:, 0], label="Susceptible", color="orange")
    ax1.plot(t, solA[:, 1], label="Infected", color="red")
    ax1.plot(t, solA[:, 2], label="Recovered", color="green")
    #plot official data
    ax1.plot(t, officialTotalInfectionsMinusOffsetScaled(), label="Official infections", color="blue")

    ax1.set_ylim([0.0, 1.0])
    ax1.set_xlim([0.0, ndays])
    ax1.set_ylabel("Ratio of Population")
    ax1.legend(loc="upper left"); ax1.grid();
    # Model B
    plt.setp(ax1.get_xticklabels(), visible=False)  # dont show xlabels on ax1
    ax2.set_xticklabels(['Sep', 'Okt', 'Nov', 'Dec', 'Jan'])

    ax2.plot(t, solB[:, 0], label="Susceptible", color="orange")
    ax2.plot(t, solB[:, 1], label="Infected", color="red")
    ax2.plot(t, solB[:, 2], label="Recovered", color="green")
    #plot official data
    ax2.plot(t, officialTotalInfectionsMinusOffsetScaled(), label="Official infections", color="blue")

    ax2.set_ylim([0.0, 1.0])
    ax2.set_xlim([0.0, ndays])
    ax2.set_ylabel("Ratio of Population")
    ax2.set_xlabel("Time")
    ax2.legend(loc="upper left"); ax2.grid();
    plt.subplots_adjust(wspace=0.1, hspace=0.09)
    plt.show()
    

def helperMap(t):
    """ Maps an input function onto another range. Used for the Half-Gaussian distribution """
    # https://math.stackexchange.com/questions/377169/going-from-a-value-inside-1-1-to-a-value-in-another-range
    a = t_0
    b = ndays
    c = -3
    d = 0
    y = (t-a) * ( (d-c) / (b-a) ) + c
    return y


def officialTotalInfectionsMinusOffset(filename="timeline-faelle-ems_period_simple.csv"):
    """
    Return a DataFrame with numbers of total_infections minus past infection at time t_0
    data source: https://www.data.gv.at/katalog/dataset/9723b0c6-48f4-418a-b301-e717b6d98c92
    """
    df = pd.read_csv("./data/"+filename, delimiter=';', decimal=',')
    df = df['total_infections']  # select column
    offset_t0 = df[0] - officialInfections_t0  # number of past infections
    df = (df.to_numpy() - offset_t0) * (1-gamma)  # scale data: subtract past infections
    return df


def officialTotalInfectionsMinusOffsetScaled(filename="timeline-faelle-ems_period_simple.csv"):
    """
    Return a DataFrame with population ratios of total_infections minus past infection at time t_0 scaled by the population N
    data source: https://www.data.gv.at/katalog/dataset/9723b0c6-48f4-418a-b301-e717b6d98c92
    """
    df = pd.read_csv("./data/"+filename, delimiter=';', decimal=',')
    # Get number of active cases (infected people) per day
    df = df['total_infections']  # select column
    offset_t0 = df[1] - officialInfections_t0  # number of past infections
    # scale data: subtract past infections, divide by population size and reduce by recovery rate
    df = (df.to_numpy() - offset_t0) / N * (1-gamma)
    return df


def plotBothSirAndOfficialTotalInfections(solutionA, solutionB, allPhaseLines):
    # Plot the 3 phase lines for S, I and R
    
    # MODEL A
    plt.figure(figsize=[8, 4])
    plt.plot(t, solutionA[:, 0], label="A:Susceptible", color="orange", marker="|")
    plt.plot(t, solutionA[:, 1], label="A:Infected", color="red", marker="|")
    if allPhaseLines: 
        plt.plot(t, solutionA[:, 2], label="A:Recovered", color="green", marker="|")

    # MODEL B
    plt.plot(t, solutionB[:, 0], label="B:Susceptible", color="orange")
    plt.plot(t, solutionB[:, 1], label="B:Infected", color="red")
    if allPhaseLines: 
        plt.plot(t, solutionB[:, 2], label="B:Recovered", color="green")

    # Plot official data
    plt.plot(t, officialTotalInfectionsMinusOffsetScaled(), label="True infected", color="blue", linewidth=3)

    plt.xticks(ticks=np.arange(0,ndays,30), labels=['Sep', 'Okt', 'Nov', 'Dec', 'Jan'])
    plt.legend(); plt.grid(); plt.title("Both SIR-models and official infection data")
    plt.xlabel("Time"); plt.ylabel("Ratio of Population")
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
    N = 9e6                         # population size
    officialInfections_t0 = 19043   # official number of recorded infected people at t0
    acc = 3                         # accuracy, number of decimal places in results

    I_0 = officialInfections_t0/N   # fraction of Infected (1353 reported cases)

    R_0 = 0.58 + 0.30        # 58% of fully Vaccinated + 30% avoiding infections

    S_0 = 1 - I_0 - R_0      # fraction of Susceptible at t0

    beta = 1.07      # transmission rate

    gamma = 1/10     # recovery rate  10 days -> 0.1

    # An infection protection of 80% is assumed (assumption copied from policy brief)
    # This level of immunity is reached through vaccination or convalescence after an infection with the virus 
    immunity_level = 0.8
    
    ##### Model A considers no decay in immunity_level over time #####
    epsilonA = 0.0

    ##### Model B considers the decay in immunity_level as a function over time #####
    # damps the decay over time
    damper = 0.01      # for 'total_infections' data

    ### Constant rate of decay
    epsilonB_constant = lambda t: 1 / ndays * t * damper

    ### sigmoid function with a damper (sigmoid is similar to a jump function)
    #epsilonB_sigmoid = lambda t: (1 / ( 1 + np.exp(-t) )) * damper

    ### Half Normal Distributed values for decay rates
    epsilonB_halfGauss = lambda t: 1. / ( np.sqrt(1 ** np.pi) ) * np.exp( -1 * np.power(helperMap(t), 2.) ) * damper

    ######### CHOICE #########
    epsilonB = epsilonB_constant


    # Print a report of the used parameters
    print("\n--------------------[ PARAMETERS ]------------------------\n")
    print(f"t_0 = {t_0}")
    print(f"timespan = {ndays} days / {ndays / 30} months")
    print(f"timestep dt = {dt} day")
    print(f"N = {int(N)}")
    print(f"S_0 = {round(S_0*100, 2)}%")
    print(f"I_0 = {round(I_0*100, 2)}%")
    print(f"R_0 = {round(R_0*100, 2)}%")
    print(f"beta = {round(beta, acc)} (transmission rate)")
    print(f"gamma = {round(gamma, acc)} (recovery rate)")
    print(f"epsilonA = {round(epsilonA, acc)} (vaccination decay rate)")
    print(f"damper = {damper} (for epsilon)")
    print(f"max reachable level of immunity = {immunity_level*100}%")


    def helperPlotBothEpsilonB(eps, xstart, xstop, damper):
        g = np.linspace(xstart, xstop, 50)
        colors=['black', 'blue']
        labels=['Constant','Half Gaussian']
        for i in range(len(eps)):
            plt.plot(g, (eps[i](g)) * 1 / damper, color=colors[i], label=labels[i])    # without damper
            #plt.plot(g, (eps[i](g)) * 1, color=colors[i], label=labels[i])     # with damper
        plt.title("Plot epsilon function for Model B\n(rate of vaccination protection decay)")
        plt.xticks(ticks=np.arange(xstart, xstop, 30), labels=['Sep', 'Okt', 'Nov', 'Dec', 'Jan'])
        #plt.yticks(np.linspace(0, 1, 11));
        #plt.ylim(0, 1);
        plt.xlim(xstart, xstop)
        plt.legend(loc='upper left')
        plt.grid(); plt.show()


    # Solve coupled system of ODEs using RK4
    solutionA = np.array(solve_ivp(
        fun=lambda t, y: SIR_modelA(y, t, beta, gamma, epsilonA, immunity_level),  # system of ODEs
        y0=[S_0, I_0, R_0],     # initial condition
        t_span=[t_0, ndays],    # time points
        t_eval=t,               # when to store computed solutions
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix

    solutionB = np.array(solve_ivp(
        fun=lambda t, y: SIR_modelB(y, t, beta, gamma, epsilonB, immunity_level),  # system of ODEs
        y0=[S_0, I_0, R_0],     # initial condition
        t_span=[t_0, ndays],    # time points
        t_eval=t,               # times when to store computed solution
        method='RK45',          # default
    ).y.T)  # select transformed solution matrix


    # Print a report of the simulation
    print("\n----------------------[ RESULT ]--------------------------\n")
    print("Day 0: 1st September, 2021")
    print(f"Day {ndays}: 31st January, 2021\n")

    pMaxInfectedModelA = round(max(solutionA[:, 1])*100, 1)
    totalMaxInfectedModelA = int(N*pMaxInfectedModelA/100)
    print(f"Infected maximum in Model A: {pMaxInfectedModelA}% / {totalMaxInfectedModelA} infected")

    pMaxInfectedModelB = round(max(solutionB[:, 1])*100, 1)
    totalMaxInfectedModelB = int(N*pMaxInfectedModelB/100)
    print(f"Infected maximum in Model B: {pMaxInfectedModelB}% / {totalMaxInfectedModelB} infected")

    pMaxInfectedReality = round(max(officialTotalInfectionsMinusOffset())/N*100,2)
    totalMaxInfectedReality = int(max(officialTotalInfectionsMinusOffset()))
    print(f"Officially infected maximum: {pMaxInfectedReality}% / {totalMaxInfectedReality} infected\n")

    percentErrorModelA = round((totalMaxInfectedReality - totalMaxInfectedModelA) / totalMaxInfectedReality * 100, 1)
    print(f"Percent error of Model A in maximum: {percentErrorModelA}%")

    percentErrorModelB = round((totalMaxInfectedReality - totalMaxInfectedModelB) / totalMaxInfectedReality * 100, 1)
    print(f"Percent error of Model B in maximum: {percentErrorModelB}%")

    print("\n[INFO] Plot is ready!")


    ################################## Plot results ####################################

    # Plot both
    helperPlotBothEpsilonB((epsilonB_constant, epsilonB_halfGauss), 0, ndays, damper)

    # Plot both models with one selected epsilon
    plotBothSirAndOfficialTotalInfections(solutionA, solutionB, False)

    # Plot both models with one sdlected epsilon
    subplotsSIR(solutionA, solutionB)
    
    # Plot both SIR models with all phase lines and the official infection data
    plotBothSirAndOfficialTotalInfections(solutionA, solutionB, True)





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                    ADDITIONAL INFORMATION 

‘RK45’ (default): Explicit Runge-Kutta method of order 5(4) [1]. 
The error is controlled assuming accuracy of the fourth-order method, 
but steps are taken using the fifth-order accurate formula (local extrapolation is done). 
A quartic interpolation polynomial is used for the dense output [2]. 
Can be applied in the complex domain.
 - Source: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp @method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""