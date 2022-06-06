# The effect of time dependent vaccination protection (SIR-model)

## Question:
How will the results of two SIR models for Covid-19 infection dynamics differ if additionally time dependence of vaccination protection is considered?  

## Hypothesis:
The model considering the time dependent vaccination protection will be more accurate. <br>

---
### Basic model: SIR (a simple epidemics model)
* N = S + I + R (the total number of individuals is constant)
* Individuals only jump between the 3 compartments. 
* Solver for system of ODEs - RK4 is very stable and accurate
---
### Assumptions: 
* Model A considers 5 months of constant vaccination protection / natural immunity	
* Model B considers the vaccination protection / natural immunity as a function of time.
* Other Assumptions:
  - **Dev Phase 1:**
    - R = Removed / Recovered / Vaccinated 
    - S -> I (Susceptiles turn infected over time)
    - I -> R (Infected recover after 10 days and gain immunity)	
  - **Dev Phase 2:**  
    - immunity duration? - 5 months 
    - R -> S (fallback after 5 months) -> How is the fallback time distributed? -> Half-Gaussian distribution
    - Birth rate b = mortality m OR both are considered negligibly small

---
