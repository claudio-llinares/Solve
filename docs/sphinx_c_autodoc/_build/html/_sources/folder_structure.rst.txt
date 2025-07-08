Folder structure
================
Here is the folder structure of the code:

| .
| ├── bin
| ├── gravity
| │   ├── commons_multigrid
| │   ├── mond_mass
| │   └── newton
| │             └── fft_solver
| ├── grid
| ├── libs
| ├── pm
| ├── tests
| └── time_evolution

Details:

- ``bin`` contains the Makefile
- ``gravity`` contains the gravity solvers.
- ``gravity/commons_multigrid`` suport routines for multigrid solvers which do not depend on the model.
- ``gravity/mond_mass`` constains the multigrid solver for the model we are dealing with (MOND with mass terms).
- ``gravity/newton`` contains an FFT based solver for Poisson's equation (needed for standard simulations).
- ``grid`` all routines related to the grid including allocation routines and initialization of the FFTW variables.
- ``libs`` miscelaneous including cosmology, utility files for allocating 3D arrays and routines for reading parameter files.
- ``pm`` (particle mesh).  Everything related to the particles at their connection to the grid.  Here we allocate memory for particles and transfer information from the particles trough the grid (for calculatign density in the grid) and viceversa (for calculating force in the particles).
- ``tests`` has routines for testing.  Here you can find the routines for the spherical test.
- ``time_evolution`` contains the main time loop.
