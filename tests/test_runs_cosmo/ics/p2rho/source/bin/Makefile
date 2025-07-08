SYSTEM="laptop"
#SYSTEM="SP6"
#SYSTEM="nekkar"

ifeq ($(SYSTEM), "laptop")
	FFTW_PREF=/home/claudio/applications/fftw-3.3.lib.openmp
        GSL_PREF=/home/claudio/applications/gsl-1.15.lib
        CC = gcc
        C = -Wall -O3 -fopenmp
endif

DEFINEFLAGS = -DOPENMP
DFLAGS = $(DEFINEFLAGS)
CFLAGS = ${C} $(DFLAGS)  -I${FFTW_PREF}/include -L${FFTW_PREF}/lib -I${GSL_PREF}/include -L${GSL_PREF}/lib
LIB    = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm
########################################################################
VPATH = ..:../libs
########################################################################
# Objects
LIBSOBJ = allocator.o cuad_min.o D_lib.o releaser.o ran3.o gasdev.o p.o p2rho_lib.o sigma8_lib.o z2t_lib.o power_lib_fftw_corrected.o 

# All objects
ALLOBJ = $(LIBSOBJ)

p2rho:   $(ALLOBJ) p2rho.o
	$(CC) $(CFLAGS) $(ALLOBJ) p2rho.o -o p2rho $(LIB)

%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@
clean:
	rm -r *.o; rm -r *~; rm p2rho
