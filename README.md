# Cloud Survival in Galactic Halo

This project is a reproduction and expansion of the original case setup outlined in - https://doi.org/10.1093/mnras/stag1405. We have successfully reproduced the original 3D setup in AGILE and included a 2D translation of the simulation. Currently, only the adiabatic case is supported.

To prescribe the background profiles for the case, we have made the use of lookup tables defined in the `mod_background` files. These must be included while building AMRVAC and, hence, custom makefiles that do so have been provided - please make use of them.

The current resolution for the lookup table is 50000 points for gravity, pressure and density. `p_offset` is a constant offset value added to the pressure profile which is a normalization choice - its value is given in CGS - $5.4 \times 10^{-13} \mathrm{dyn/cm}^2$ and has been confirmed by the authors of the paper.

While running the simulation please adjust the maximum number of blocks based on your compute. 

## Explanation of the refinement routine

We make use of the native tracer support in AMRVAC that advects a scalar field along with the solving for the other variables. This is later used to track cloud particles and actively used to perform refinement during the simulation.

In the `usr_amr_refine` subroutine, we first identify whether the block contains cloud material exceeding a certain threshold value of the tracer. The second check also identifies whether enough density of the cloud material is present to prevent unnecessary over-refinement of every block containing the tracer. 

These threshold values are as follows
1. `C_min` - Currently this is taken as `1d-3` as testing with even finer options gives the same result.
2. `rho_bulk` - Identifies the bulk of the cloud material to refine - currently at $0.01$ times the original cloud density. Can. be tested with coarser values such as $0.1$ which give the same results

```
subroutine usr_amr_refine(igrid, level, ixI^L, ixO^L, qt, w, x, &
                            refine, coarsen)
    .
    .
    .

    has_tracer = any(w(ixO^S, tracer(1)) > C_min * w(ixO^S, rho_))
    has_bulk   = any(w(ixO^S, tracer(1)) > C_min * w(ixO^S, rho_) .and. &
                     w(ixO^S, rho_)      > rho_cl_bulk)

    if (has_bulk) then
       refine  = 1                           ! bulk cloud: force refinement
       coarsen = 0
    else if (has_tracer) then
       refine = -1                           ! diluted remnant: stop refining it
       if (level > level_bg) coarsen = 1     ! and release the cells
    end if
    ! no tracer at all -> refine = coarsen = 0 -> Lohner decides, untouched
  end subroutine usr_amr_refine
```

Since, the background can also be sensitive to refinement, we have defined an upper threshold to not refine it beyond that. Currently this is at $2$ but testing has showed that increasing it does not change the results.
