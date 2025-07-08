reset
unset multiplot

set autoscale

set size square
set grid


plot "../run_0064_0064mpc_01_std/conv.0.0" \
     binary format="%7float64" u \
     ($1-32.0):($2-32.0) w d lc 1
