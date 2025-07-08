SYSTEM="laptop"

DEFINEFLAGS = -DOPENMP
DFLAGS = $(DEFINEFLAGS)

ifeq ($(SYSTEM), "laptop")
	FFTW_PREF=/home/claudio/applications/fftw-3.3.lib.openmp
	GSL_PREF=/home/claudio/applications/gsl-1.15.lib
	CC = gcc
	C = -Wall -fopenmp -O3
	LIB = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm
	CFLAGS = ${C} $(DFLAGS) -I${FFTW_PREF}/include -L${FFTW_PREF}/lib -I${GSL_PREF}/include -L${GSL_PREF}/lib
endif


########################################################################
VPATH = ..:../gravity:../gravity/newton:../libs:../grid:../pm
########################################################################
# Objects
OBJS = init_physics.o

LIBSOBJ =  utils.o read_params.o allocator.o releaser.o D_lib.o D_dot_lib.o cosmology.o

GRIDOBJ = init_grid.o init_fftw.o clean_up_grid.o

PARTOBJ = init_parts.o kernels.o interp_force.o write_particles.o zeldovich.o

GRAVITYOBJ = solve_poisson_fftw.o get_forces.o calc_newt_source.o write_initial_potential.o

# All objects
ALLOBJ = $(OBJS) $(LIBSOBJ) $(GRIDOBJ) $(PARTOBJ) $(GRAVITYOBJ) $(TIMEEV)

solve:   $(ALLOBJ) solve.o
	$(CC) $(CFLAGS) $(ALLOBJ) solve.o -o solve $(LIB)

%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@ #-pg
clean:
		rm *.o; cd ../; find . -type f -name "*~" -delete; cd bin
