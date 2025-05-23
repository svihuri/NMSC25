## Stable Fluid Solver

Final project for NMSC 2025 by Saana Vihuri. Modular C++ version of Jos Stam's stable fluid solver https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf. <br>

Simulation sample, Kelvin-Helmholtz instability in 2D:

![](kevin_helmholz_instability.gif)

### How to compile and run: <br>

git clone https://github.com/svihuri/NMSC25.git <br>
cd NMSC25/ <br>

 g++ main.cpp -O3 -lfftw3f -o main <br>
 ./main <br>
 python3 plot_simulation.py

 OR

 ./reproduce_simulation.sh

#### Dependencies 

Tested on: <br>

Ubuntu 24.04 <br>
g++ 13.3.0 <br>
FFTW 3.3.10 <br>

No fftw3.h --> sudo apt-get install libfftw3-dev 

Plotting:  <br>

Python 3.12.3 <br>
pandas 2.1.4 <br>
numpy 1.26.4 <br>

ModuleNotFoundError: No module named 'pandas' --> sudo apt-get install python3-pandas or use pip <br>
