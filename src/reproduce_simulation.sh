#! /bin/bash

g++ main.cpp -O3 -lfftw3f -o main
./main
python3 make_animation.py
python3 plot_snapshots.py
