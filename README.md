Supplementary Material for "Ecological inheritance facilitates the coexistence of environmental helpers and free-riders"


The files "sims_island_fixed_patch_size.jl", "sims_island_stoch_patch_size.jl", "sims_lattice_fixed_patch_size.jl" contain the simulation programs used in the study to simulate the evolution of the investment z into a durable environment. The simulations are written using the programming language julia (version 1.9.2). With julia installed, the programs can then be launched from the working directory using the command line: julia program_name. Each simulation program is annotated.

-the program "sims_island_fixed_patch_size.jl" simulates evolution in a population divided into a very large number of patches uniformly connected by dispersal  (that is, under the island model of dispersal), where each patch has a fixed size N. The procedure is described in Appendix C-7.
-the program "sims_island_stoch_patch_size.jl" simulates evolution in a population divided into a very large number of patches uniformly connected by dispersal, where the number of individuals per patch fluctuates due to density-dependent survival. The procedure is described in Appendix D.
-the program "sims_lattice_fixed_patch_size.jl" simulates evolution in a population divided into a very large number of patches distributed in a two-dimensional lattice, connected by stepping-stone dispersal (that is, offspring can only disperse to neighbouring patches). Each patch has a fixed size N. The procedure is described in Appendix E-2.


The files "Supplementary_Movie_m_0.8_lambda_0.3.gif" and "Supplementary_Movie_m_0.8_lambda_0.9.gif" contain the polymorphism-driven spatio-temporal dynamics of the environment, each snapshot representing the environmental quality of all patches in the population. These follow the spatio-temporal dynamics of the environment for 100 generations. The simulations used to represent these dynamics are the same ones used to draw Figure 5 of the main text (see the legend of Figure 5 for parameter values). 
