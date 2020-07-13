# Hybrid-RSP

The matlab code in this repository is the code used for the simulations for the examples in the paper "A Hybrid Control Algorithm for Gradient-Free Optimization using Conjugate Directions" by Alessandro Melis, Ricardo G. Sanfelice, and Lorenzo Marconi.

It implements a robust hybrid controller based on conjugate directions able to minimize a rather general objective function.

To run the code, the Hybrid Systems Toolbox is needed. It can be downloaded at https://www.mathworks.com/matlabcentral/fileexchange/41372-hybrid-equations-toolbox-v2-04.

In order to run the simlations, simply run the file "Run_this_file".

The file "initialization" defines all the initial conditions of the states, as well as design parameters, for the simulation. They can both be directly defined, or randomly generated.
In this file, 3 modes of operation have been considered:
- Mode 1, considering no noise acting on the measurements of the objective function;
- Mode 2, considering noise acting on the measurements of the objective function. In particular, the noise generated is always able to steer the optimization variable X away from the set of minima of the objective function;
- Mode 3, considering noise acting on the measurements of the objective function. In this case, by imposing a lower bound on the global step size Phi, robust convergence to a neighborhood to the set of minima.

For any information or problem, please write at alessandro.melis4@unibo.it.

Best Regards,
Alessandro Melis, Ricardo G. Sanfelice, Lorenzo Marconi
