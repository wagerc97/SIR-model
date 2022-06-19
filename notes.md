# Notes 
IMPORTANT: 
* [Model Structure.xlsx](C:\Users\wager\OneDrive\UNI_WIEN\1LV_online\4thSem\PDE\Model Structure.xlsx) for more information and sources on the parameters
* [PDE Roadmap.docx](C:\Users\wager\OneDrive\UNI_WIEN\1LV_online\4thSem\PDE\PDE Roadmap.docx)   

## Missing model parameters
- ``R_t0`` the recovered fraction of the population 
  (``S_t0`` would follow from ``R_t0``)
- ``beta / a / R(eff)`` the infection rate (effective reproduction rate)
  - ``beta`` = 1.07 according to [https://orf.at/corona/daten/oesterreich](https://orf.at/corona/daten/oesterreich)
  - ``beta`` could also be defined by ``beta = k * q * D `` (k = contact rate, q = probability of an infection, D = 7 duration of infectious state in days)
    But we don't have accurate sources for that.
- ``damper`` = 0.05 to damp the effect of the ``epsilonB`` which models the decay of the vaccination  
- ``R0 = beta/gamma`` so infection over recovery rate

## Offene Fragen
- _Wie passe ich die Zahlen an?_ **damper angepasst **
- Skaliert auf 9Mio Einwohner in Ö sind das SEHR kleine Zahlen...?
- _Soll ich R(eff) die effektive Reproduktionsrate plotten? R(eff) wäre ein schönes Vergleichsmaß zu realen Daten!_ **NEIN, zu viel Aufwand**
- _Soll ich S_crit über die Zeit plotten?_ **NEIN, zu viel Aufwand**

## TODO
- _Ich muss die realen Daten aus dieser Zeit finden! -> Plot_ **FERTIG**

## Release Notes
10.Juni 2022: Plakatives Beispiel mit Plots die zeigen, wie sich die Annahme auf ein solches System auswirtk   
19.Juni 2022: Initial values und damper angepasst an echte Daten aus offiziellen Quellen.  
- Resultat sind Plots ohne Recovered population ratio
- eine Kurve, die sehr nahe an die echten Daten herankommt: 
  - gleiche form 




