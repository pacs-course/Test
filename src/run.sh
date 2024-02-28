#!/bin/bash
make
./main_testLMat
gnuplot -p <<EOF 
plot "result.dat" w lp title "uh"
EOF
