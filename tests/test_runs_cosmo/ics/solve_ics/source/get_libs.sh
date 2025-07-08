#!/bin/sh

cp ../../../varios/allocator/allocator_mio/allocator.c libs
cp ../../../varios/allocator/allocator_mio/allocator.h libs

cp /home/claudio/source/numerics/recipies/recipies_original/src/gasdev.c libs 
cp /home/claudio/source/numerics/recipies/recipies_original/src/ran1.c libs 

# For the analytic densities and Zeldovich
cp ../../../cosmo_util/D/D_lib.c libs
cp ../../../cosmo_util/D/D_lib.h libs
cp ../../../cosmo_util/D/D_dot_lib.c libs
cp ../../../cosmo_util/D/D_dot_lib.h libs

# For Bob's power spectra
# WARNING! THIS ROUTINES WERE MODIFIED IN THIS DIRECTORY.  
# DO NOT COPY THEM AGAIN.
##cp ../../p2rho/p2rho_lib.c libs
##cp ../../p2rho/p2rho_lib.h libs

