Data structure
==============

We have the following structures defined in the code:

- ``params``: defined in solve.h.  It contains input parameters and some physical constants (including expansion factor).
- ``parts``: defined in ``pm/init_parts.h``.  It contains particles' positions and velocities, which are allocated in ``pm/init_parts.c``
- ``grids``: defined in ``grid/init_grid.h``.  If contains all the fields and space for doing FFTs.  Eveything is allocated in ``grid/init_grid.h``.  Note that there are additional grids allocated by the multigrid solver.
- ``grids_gravity_mm``: defined in ``gravity/mond_mass/solve_multigrid_mm.h``.  It contains everything related to the MOND-MM multigrid solver (including the coarse grids, which are allocated in ``gravity/mond_mass/solve_multigrid_mm.c``.


Detailed definitions follow:

.. autoctype:: solve.h::params
.. autoctype:: pm/init_parts.h::particles
.. autoctype:: grid/init_grid.h::grids
.. autoctype:: gravity/mond_mass/solve_multigrid_mm.h::grids_gravity_mm



