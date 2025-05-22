Final project for NMSC 2025 by Saana Vihuri. Modular C++ version of Jos Stam's stable fluid solver https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf. <br>

Simulation sample: Kelvin-Helmholtz instability.

![](kevin_helmholz_instability.gif)

How to compile and run: <br>
 
 g++ main.cpp -O3 -lfftw3f -o main <br>
 ./main <br>
 python3 plot_simulation.py
