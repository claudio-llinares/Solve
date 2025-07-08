reset
unset multiplot

set size square
set grid

set terminal x11 
set multiplot layout 2,2

plot "../run_0064_0064mpc_newt/conv.0.0" \
     binary format="%7float64" u \
     ($1-32.0):($2-32.0) w d lc 1

plot "../run_0064_0064mpc_mond_small_a0_without_m/conv.0.0" \
     binary format="%7float64" u \
     ($1-32.0):($2-32.0) w d lc 1

plot "../run_0064_0064mpc_mond_without_m/conv.0.0" \
     binary format="%7float64" u \
     ($1-32.0):($2-32.0) w d lc 1

plot "../run_0064_0064mpc_mond_with_m/conv.0.0" \
     binary format="%7float64" u \
     ($1-32.0):($2-32.0) w d lc 1
