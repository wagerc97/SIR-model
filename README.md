## The effect of time dependent vaccination protection (SIR-model)

University of Vienna, S2022  
250070 SE Seminar Applied PDE (2022S)  
Supervisor: Univ.-Prof. Dr. Norbert Mauser  
View the **[Report](./report_SIR-model_WAGER.pdf)** to read the whole study! 


### Question:
How will the results of two simple SIR models for Covid-19 infection dynamics differ if additionally time dependence of vaccination protection is considered?  

### Hypothesis:
The model considering the time dependent vaccination protection will be more accurate. <br>

### Model simulation: 
#### Basic model: SIR (a simple epidemics model)
* _N_ is a homogenously mixed population with size _N_=1   
  _N = S + I + R_
* Individuals only jump between the 3 compartments
* RK4 solver used for ODEs is very stable and accurate

#### The two models: 
* Model A considers 5 months of constant vaccination protection / natural immunity.	
* Model B considers the vaccination protection / natural immunity as a function of time.
* A time frame of September 1st, 2021 until January 31st, 2021 is studied.

#### Assumptions: 
- The average convalescence period of a Covid-19 infection is 10 days
- The maximum level of immunity _e_ = 80% on average
- There are no new vaccinations during the studied time period
- A damping factor _f_ = 0.01 is used to calibrate the effect of the central model B assumption
- Constant transmission rate _beta_ =1.07 
- The population ratio of recovered (= immune) is calibrated to be 0.88 
by adding an estimated amount of 0.3

#### Acknowledgement
This is version 2.0 of the SIR-model and report.  
I thank Professor Norbert Mauser for the greate support.  
The final version incorporates recommendations/critique from Professor Dieter Suess. I thank you for the input.
