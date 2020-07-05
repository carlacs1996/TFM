# Master's project code

This repo contains all code needed to reproduce the experiments 
shown in my Master's project.

potential.m 
    The code solves ODE (3.3) using Euler method and plots 
    the position of the particle over time (Figure 3.2).

periodic.m
    The code solves ODE (3.8) using Euler method and plots 
    the position of the particle over time (Figure 3.3). 

white.m
    The code solves ODE (3.28) using Euler Maruyama method and plots 
    the position of the particle over time for three different values 
    of kappa and computes the mean exit time for each case (Figure 3.4). 

finalSDE.m
    The code solves ODE (3.34) using Euler Maruyama method and plots 
    the position of the particle over time for three different values 
    of kappa (Figure 3.5). 

SR1.m 
    The code plots the position of the particle and its periodogram 
    for a single realisation and an an average for both the (only) 
    noise SDE (3.28) and our final SDE (3.34) (Figure 3.6).

SR2.m
    The code plots the periodogram for different values of the weak 
    force frequency omega_s (Figure 3.7).

SR3.m
    The code plots the position of the particle (a single realisation 
    and an average) and its periodogram for different values of the noise 
    strength kappa (Figure 3.8).