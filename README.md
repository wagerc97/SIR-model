## SIR  - an Epidemic Model
**Question**:How will the results of two SIR models for Covid-19 infection dynamics differ if one additionally considers time dependence of vaccination protection?  <br>
**Hypothesis:** The model considering the time dependent vaccination protection will be more accurate.  <br>
---
#### Basic model: SIR (a simple epidemics model)
* N = S + I + R (the total number of individuals is constant)
* Individuals only jump between the 3 compartments. 
* Solver for system of ODEs - RK4 is very stable and accurat
---
#### Assumptions: 
* Model A considers 5 months of constant vaccination protection / natural immunity	
* Model B considers the vaccination protection / natural immunity as a function of time over the course of 5 months.
* Other Assumptions:
  - R = Removed / Recovered / Vaccinated 
  - S -> I (Susceptiles turn infected over time)
  - I -> R (Infected recover after 10 days and gain immunity)	Phase2: immunity duration? - 3 / 5 / 6 months 
  - R -> S (fallback after 5 months)
  - Birth rate b = mortality m -> negligibly small
---
