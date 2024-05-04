#!/bin/bash

cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print($2)}' > rhoee.dat
cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print($3)}' > rhoxx.dat
cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print($6)}' > rhoeebar.dat
cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print($7)}' > rhoxxbar.dat
cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print(sqrt($4*$4+$5*$5))}' > rhoex.dat
cat mu100-400by400by100_Ev0_azym_grid.sum | awk '{print(sqrt($8*$8+$9*$9))}' > rhoexbar.dat
cat temporary.raw | head -n 1 > test_v0_init.dat
cat temporary.raw | tail -n 1 > test_v0.dat

cat temporary.raw | head -n 120 | tail -n 1 > test_v0_120.dat
cat temporary.raw | head -n 180 | tail -n 1 > test_v0_180.dat
cat temporary.raw | head -n 200 | tail -n 1 > test_v0_200.dat
cat temporary.raw | head -n 300 | tail -n 1 > test_v0_300.dat
cat temporary.raw | head -n 400 | tail -n 1 > test_v0_400.dat
cat temporary.raw | head -n 410 | tail -n 1 > test_v0_410.dat
cat temporary.raw | head -n 500 | tail -n 1 > test_v0_500.dat
cat temporary.raw | head -n 290 | tail -n 1 > test_v0_nXXX.dat
cat temporary.raw | head -n 600 | tail -n 1 > test_v0_600.dat
cat temporary.raw | head -n 700 | tail -n 1 > test_v0_700.dat
cat temporary.raw | head -n 800 | tail -n 1 > test_v0_800.dat
cat temporary.raw | head -n 900 | tail -n 1 > test_v0_900.dat

cat mnr.sum | awk '{print($1)}' > time.dat
cat mnr.sum | awk '{print($2)}' > mnr_nu.dat
cat mnr.sum | awk '{print($3)}' > mnr_nu0.dat
cat mnr.sum | awk '{print($4)}' > mnr_lam.dat
cat mnr.sum | awk '{print($5)}' > mnr_vac.dat

##### off-diagonal ex ell components (norm)
##cat moments.sum | awk '{print($1)}' > time.dat
#cat moments.sum | awk '{print($2)}' > N0.dat
#cat moments.sum | awk '{print($3)}' > N1.dat
#cat moments.sum | awk '{print($4)}' > N2.dat
#cat moments.sum | awk '{print($5)}' > N3.dat
#cat moments.sum | awk '{print($6)}' > N4.dat
#cat moments.sum | awk '{print($7)}' > N5.dat

#### diagonal xx ell components
#cat moments.sum | awk '{print($8)}' > N0xx.dat
#cat moments.sum | awk '{print($9)}' > N1xx.dat
#cat moments.sum | awk '{print($10)}' > N2xx.dat
#cat moments.sum | awk '{print($11)}' > N3xx.dat
#cat moments.sum | awk '{print($12)}' > N4xx.dat
#cat moments.sum | awk '{print($13)}' > N5xx.dat

