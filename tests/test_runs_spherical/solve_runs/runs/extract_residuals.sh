#!/bin/bash


grep verbose_multi run_032_016mpc/run.log | cut -f2 -d= > residuals.dat


