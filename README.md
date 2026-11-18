# Cloud Survival in Galactic Halo

This project is a reproduction and expansion of the original case setup outlined in - https://doi.org/10.1093/mnras/stag1405. We have successfully reproduced the original 3D setup in AGILE and included a 2D translation of the simulation. Currently, only the adiabatic case is supported.

To prescribe the background profiles for the case, we have made the use of lookup tables defined in the `mod_background` files. These must be included while building amrvac and, hence, custom makefiles that do so have been provided - please make use of them.

The current resolution for the lookup table is 4000 points for gravity, pressure and density. `p_offset` is a constant offset value added to the pressure profile which is a normalization choice - its value is given in CGD dyn/cm^2 and has been confirmed by the authors of the paper.

While running the simulation please adjust the maximum number of blocks based on your compute. 