#!/bin/bash

SOURCE_DIR=../../../source
SOURCE_LOCAL=./source

mkdir $SOURCE_LOCAL
mkdir $SOURCE_LOCAL/bin
mkdir $SOURCE_LOCAL/gravity
mkdir $SOURCE_LOCAL/gravity/mond_mass
mkdir $SOURCE_LOCAL/gravity/newton
mkdir $SOURCE_LOCAL/gravity/newton/fft_solver
mkdir $SOURCE_LOCAL/grid
mkdir $SOURCE_LOCAL/libs
mkdir $SOURCE_LOCAL/pm
mkdir $SOURCE_LOCAL/tests
mkdir $SOURCE_LOCAL/time_evolution

cp $SOURCE_DIR/bin/Makefile $SOURCE_LOCAL/bin
cp $SOURCE_DIR/gravity/commons_multigrid/prolonge.c $SOURCE_LOCAL/gravity
cp $SOURCE_DIR/gravity/commons_multigrid/restrict.c $SOURCE_LOCAL/gravity
cp $SOURCE_DIR/gravity/get_forces.c $SOURCE_LOCAL/gravity
cp $SOURCE_DIR/gravity/get_gravity.c $SOURCE_LOCAL/gravity
cp $SOURCE_DIR/gravity/get_lap.c $SOURCE_LOCAL/gravity
cp $SOURCE_DIR/gravity/mond_mass/cycles_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/differential_operators_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/go_up_and_down_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/gs_residual_truncation_error_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/initial_guess_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/solve_multigrid_mm.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/solve_multigrid_mm.h $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/mond_mass/solve_spherical_mond.c $SOURCE_LOCAL/gravity/mond_mass
cp $SOURCE_DIR/gravity/newton/calc_newt_source.c $SOURCE_LOCAL/gravity/newton
cp $SOURCE_DIR/gravity/newton/fft_solver/solve_poisson_fftw.c $SOURCE_LOCAL/gravity/newton/fft_solver
cp $SOURCE_DIR/grid/add_two_fields.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/clean_up_grid.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/get_mean_max_two_fields.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/init_fftw.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/init_grid.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/init_grid.h $SOURCE_LOCAL/grid
cp $SOURCE_DIR/grid/write_phi.c $SOURCE_LOCAL/grid
cp $SOURCE_DIR/init_physics.c $SOURCE_LOCAL
cp $SOURCE_DIR/init_physics.h $SOURCE_LOCAL
cp $SOURCE_DIR/libs/allocator.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/a_of_time.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/cosmology.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/ran3.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/read_params.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/releaser.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/libs/utils.c $SOURCE_LOCAL/libs
cp $SOURCE_DIR/pm/calc_delta_cic.c $SOURCE_LOCAL/pm
cp $SOURCE_DIR/pm/init_parts.c $SOURCE_LOCAL/pm
cp $SOURCE_DIR/pm/init_parts.h $SOURCE_LOCAL/pm
cp $SOURCE_DIR/pm/input_particles.c $SOURCE_LOCAL/pm
cp $SOURCE_DIR/pm/interp_force_cic.c $SOURCE_LOCAL/pm
cp $SOURCE_DIR/pm/write_particles.c $SOURCE_LOCAL/pm
cp $SOURCE_DIR/solve.c $SOURCE_LOCAL
cp $SOURCE_DIR/solve.h $SOURCE_LOCAL
cp $SOURCE_DIR/tests/static_test.c $SOURCE_LOCAL/tests
cp $SOURCE_DIR/time_evolution/time_evolution_standard_variable.c $SOURCE_LOCAL/time_evolution

sphinx-c-apidoc -f -o ./output $SOURCE_LOCAL


