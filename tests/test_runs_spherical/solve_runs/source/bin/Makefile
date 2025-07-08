SYSTEM="laptop"
#SYSTEM="cosma"

DEFINEFLAGS = -DOPENMP
DFLAGS = $(DEFINEFLAGS)

ifeq ($(SYSTEM), "laptop")
	FFTW_PREF=/home/claudio/applications/fftw-3.3.lib.openmp
	GSL_PREF=/home/claudio/applications/gsl-1.15.lib
	CC = gcc
#	C = -Wall -fopenmp -std=gnu99  -g  #-Ofast #-O3 -g -pg
	C = -Wall -fopenmp -std=gnu99 -Ofast #-O3 -g -pg
	#C = -Wall -fopenmp -std=gnu99  -Ofast #-g  #-Ofast #-O3 -g -pg
	LIB = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm
	CFLAGS = ${C} $(DFLAGS) -I${FFTW_PREF}/include -L${FFTW_PREF}/lib -I${GSL_PREF}/include -L${GSL_PREF}/lib
endif
ifeq ($(SYSTEM), "cosma")
	FFTW_PREF=#/home/claudio/applications/fftw-3.3.lib.openmp
	GSL_PREF=#/home/claudio/applications/gsl-1.15.lib
	CC = icc
	C = -Wall -fopenmp -O3 -ipo #-fp-model  fast=2 -xHost -no-prec-div -static
	LIB = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm
	CFLAGS = ${C} $(DFLAGS) -I${FFTW_PREF}/include -L${FFTW_PREF}/lib -I${GSL_PREF}/include -L${GSL_PREF}/lib
endif


########################################################################
VPATH = ..:../libs:../grid:../pm:../time_evolution:../tests:../gravity:../gravity/newton:../gravity/newton/fft_solver:../gravity/commons_multigrid:../gravity/mond_mass
########################################################################
# Objects
OBJS = init_physics.o

LIBSOBJ = utils.o read_params.o allocator.o releaser.o cosmology.o a_of_time.o ran3.o

GRIDOBJ = init_grid.o init_fftw.o clean_up_grid.o write_phi.o get_mean_max_two_fields.o add_two_fields.o

PARTOBJ = init_parts.o calc_delta_cic.o interp_force_cic.o input_particles.o write_particles.o

GRAVITYOBJ = solve_poisson_fftw.o get_forces.o calc_newt_source.o get_lap.o get_gravity.o

GRAVITYMULTIGRID = prolonge.o restrict.o

GRAVITYMONDMASS = solve_multigrid_mm.o initial_guess_mm.o cycles_mm.o gs_residual_truncation_error_mm.o differential_operators_mm.o go_up_and_down_mm.o

TIMEEV = time_evolution_standard_variable.o

TESTSOBJ = static_test.o

# All objects
ALLOBJ = $(OBJS) $(LIBSOBJ) $(GRIDOBJ) $(PARTOBJ) $(GRAVITYOBJ) $(TIMEEV) $(TESTSOBJ) $(GRAVITYMULTIGRID) $(GRAVITYMONDMASS)

solve:   $(ALLOBJ) solve.o
	$(CC) $(CFLAGS) $(ALLOBJ) solve.o -o solve $(LIB)

%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@ #-pg
clean:
		rm *.o; cd ../; find . -type f -name "*~" -delete; cd bin
