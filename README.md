Work in progress version of code to run closed-loop simulations of Kukskokwim "like" system to evaluate the performance of alternative management procedures in the face of environmental stochasticity and uncertainty in drivers of recruitment:

>Connors B. M., B. Staton, L. Coggins, C. Walters, M. Jones, D. Gwinn, M. Catalano and S. Fleischman. In preparation. Incorporating harvest-population diversity trade-offs into salmon management in large river basins: a case study of Kuskokwim River Chinook.

Posterior samples required to parameterize the closed loop simulations are available upon request from brendan.connors@dfo-mpo.gc.ca and are from analyses described in:

>Staton B., M. Catalano, B.M. Connors, L. Coggins, M. Jones, C. Walters, S. Fleischman and D. Gwinn. In preparation. Evaluation of methods for quantifying population diversity in mixed-stock Pacific salmon fisheries.



## Files
- `make.R`: All code required to reproduce current set of simulations and figures.

- `load.R`: Loads packages and scripts necessary for analysis. This file should be sourced prior to running other scripts.

- `functions.R`: All functions written for the analysis should be placed in this file.
  
- `close_loop_sims.R`: Run closed loop forwad simulations.

- `figures.R`: Generate figures not associated with simulations (e.g., equilibrium trade-offs, populaiton diversity )
  
- `Simulation_summary.Rmd`: R Markdown doc that summarizes simulations, and generates figures for manuscript (in `images` folder)

- `ricker_to_bevholt.R`: a simulation of the system under Beverton-Holt and Ricker stock-recruitment formulations, and  a range of alpha multipliers for the BH formulation, to estimate what the multiplier should be to generate equivilent equilibrium conditions (i.e., escapements). This allows comparisons of the trade-offs that emerge under the Beverton-Holt and Ricker stock-recruitment formulations to not be confounded by differences in equilibrium states.


