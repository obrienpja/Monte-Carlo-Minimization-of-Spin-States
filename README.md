#Monte Carlo minimization of spin states

This code was used as part of a simulated annealing procedure to find the ground state of a classical Heisenberg model on the kagome lattice. We used perturbation theory on the kagome Kondo lattice model in order to obtain the Heisenberg model. The ground state phase diagram was found and reported in the Physical Review B publication, https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.024401.

The simulated annealing procedure we used was to run this Monte Carlo code for successively lower temperatures. In each Monte Carlo run, the parameter eta was adjusted in addition to the temperature, in order to achieve an acceptance rate of iterations of 50%.
