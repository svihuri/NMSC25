#include <cmath>
#include <complex>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <cstring> // for memset
#include <iomanip> 
#include <string>
#include <algorithm> 


constexpr float PI = 3.14159265358979323846f;

/* 
2D FluidSolver by Saana Vihuri (2025) based on Jos Stam's stable fluid solver (https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/jgt01.pdf).
The program implements Navier-Stokes equations for incompressible fluids.
Viscous diffusion and enforsing conservation of mass are done in frequency domain using FFTW3 library.

Building blocks for the solver (add forces, advect, diffuse, project) are taken from the Stam's algorithm.
Main changes:
1. Class FluidSolver is implemented to encapsulate the solver and its parameters.
2. Class constructor and destructor handle memory allocation and deallocation for FFTW.
3. Arrays and other data structures are not reused and named according to their purpose.
4. FFT is not done in-place

to do: make non-copyable and non-movable? 
note: Stam's code uses single precision floats so we use float instead of double
*/

class FluidSolver {
public:
    FluidSolver(const int gridSize, const float viscosity, const float timestep, std::ofstream& logfile)
        : n(gridSize), visc(viscosity), dt(timestep), logFile(logfile) {
        /* 
        parameters:
        n: grid size (n x n)
        visc: viscosity
        dt: time step
        logFile: log file stream
        */

        // FFT transformation is not done in-place, so no need for padding
        int size = n * n;
        int fft_size = n * (static_cast<int>(n / 2) + 1);       // fftw frequency domain output is (n x (n/2+1)) complex array, this affects some loops 

        // following Stam's code, u ~v_x, v ~v_y
        // allocate memory for velocity fields 
        u_now = (float*)fftwf_malloc(sizeof(float) * size);        // current velocity
        v_now = (float*)fftwf_malloc(sizeof(float) * size);
        u_prev = (float*)fftwf_malloc(sizeof(float) * size);       // previous step velocity
        v_prev = (float*)fftwf_malloc(sizeof(float) * size);

        // allocate buffers for FFT
        u_freq = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fft_size);    
        v_freq = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * fft_size);
        u_fft = (float*)fftwf_malloc(sizeof(float) * size);        
        v_fft = (float*)fftwf_malloc(sizeof(float) * size);

        // initialize current velocity fields to zero
        std::fill(u_now, u_now + size, 0.0f);
        std::fill(v_now, v_now + size, 0.0f);

        // initialise dye field (for visualisation)
        D_now.resize(size);
        D_prev.resize(size);

        // create FFTW plans (initialisation)
        // in and out arrays are specified as the transformation is not in-place
        u_plan_fwd = fftwf_plan_dft_r2c_2d(n, n, u_fft, u_freq, FFTW_ESTIMATE);
        v_plan_fwd = fftwf_plan_dft_r2c_2d(n, n, v_fft, v_freq, FFTW_ESTIMATE);
        u_plan_bwd = fftwf_plan_dft_c2r_2d(n, n, u_freq, u_fft, FFTW_ESTIMATE);
        v_plan_bwd = fftwf_plan_dft_c2r_2d(n, n, v_freq, v_fft, FFTW_ESTIMATE);

        // allocate for external forces
        f_x = (float*)fftwf_malloc(sizeof(float) * size);
        f_y = (float*)fftwf_malloc(sizeof(float) * size);
        
        // initialise to zero, default setting
        std::fill(f_x, f_x + size, 0.0f);
        std::fill(f_y, f_y + size, 0.0f);

        writeLog(logFile, "FluidSolver initialized successfully");
        
    }

    ~FluidSolver() {                        
        // destructor to clean up allocated memory & FFTW plans
        fftwf_destroy_plan(u_plan_fwd);
        fftwf_destroy_plan(v_plan_fwd);
        fftwf_destroy_plan(u_plan_bwd);
        fftwf_destroy_plan(v_plan_bwd);

        fftwf_free(u_now);  fftwf_free(v_now);
        fftwf_free(u_prev); fftwf_free(v_prev);
        fftwf_free(u_freq); fftwf_free(v_freq);
        fftwf_free(u_fft);  fftwf_free(v_fft);
        fftwf_free(f_x);    fftwf_free(f_y);
    }
     
     void solve() {
        /*solve runs one simulation step and has 4 main parts:
        1. addForces: add external forces to the velocity fields (spatial domain)
        2. advect: semi-Lagrangian advection step (spatial domain)
        3. performFFT: do FFT on velocity fields depending on direction  
        4. applyViscosityAndProject: apply viscosity and project velocity field to conserve mass (frequency domain)
        */
        if (count == 1) {
            writeLog(logFile, "Solver started");
        }
        addForces();
        advect();
        copyToFFT(u_now, u_fft);
        copyToFFT(v_now, v_fft);
        performFFT(FFTDirection::Forward, u_fft, u_freq, u_plan_fwd);
        performFFT(FFTDirection::Forward, v_fft, v_freq, v_plan_fwd);
        applyViscosityAndProject();
        performFFT(FFTDirection::Inverse, u_fft, u_freq, u_plan_bwd);
        performFFT(FFTDirection::Inverse, v_fft, v_freq, v_plan_bwd);
        normalize();
        count++;

    }
    void setExternalForces(int mode, float U0, float delta, float A, float sigma, float k) {
        // if not called, default is zero
        // mode 0: set external forces to zero, 1: set external forces from assignment
        if (mode == 0) {
            std::fill(f_x, f_x + n * n, 0.0f);
            std::fill(f_y, f_y + n * n, 0.0f);
        }
        else if (mode == 1) {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    float y_j = ((float)j+0.5f) / n;
                    float x_i = ((float)i+0.5f) / n;
                    f_x[i + n * j] = U0 * tanhf((y_j - 0.5f)/delta); 
                    f_y[i + n * j] = A * sinf(2*PI*k*x_i) * expf(-((y_j-0.5f)*(y_j-0.5f))/(2*sigma*sigma));     
            }
        }
        writeLog(logFile, "External forces set with mode " + std::to_string(mode));
        
    }
    void setInitialDyeField() {
        // set initial dye field (for visualisation)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                float y_j = ((float)j+0.5f) / n;
                if (y_j < 0.5f) {
                    D_now[i + n * j] = 0.0f; 
                } else {
                    D_now[i + n * j] = 1.0f; 
                }
            }
        }
        writeLog(logFile, "Initial dye field set");
            }

    void writeToFile(std::ofstream& file, float time, int step) {
        // write simulation output to file
        // format: time; step; i; j; u; v; D where i,j are grid coordinates
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int idx = i + n * j;
                file << std::fixed << std::setprecision(4)
                     << time << ";" << step << ";" << i << ";" << j << ";"
                     << u_now[idx] << ";" << v_now[idx] << ";" << D_now[idx] << "\n";
            }
        }
        // writing timestep
        std::ofstream timestepfile;
        timestepfile.open("timestep", std::ios::app);

        timestepfile << time << "\n";

        timestepfile.close();

        std::ofstream Dfile;
        Dfile.open("D", std::ios::app);

        for (int i = 0; i<n; i++) {
            for (int j = 0; j<n; j++) {
                Dfile << D_now[i + n * j] << " ";
            }
        }

        Dfile << "\n";
        Dfile.close();
    }

private:
    int n;                     
    float visc;                
    float dt;   
    std::ofstream& logFile;               
    int padded_size;           

    // parameters for external forces
    float* U0;
    float* delta;
    float* A;
    float* sigma;
    float* k;
    
    // real-space buffers
    float* u_now;
    float* v_now;
    float* u_prev;
    float* v_prev;
    float* u_fft;             
    float* v_fft;            
    float* f_x;            
    float* f_y;             
    
    // dye field 
    std::vector<float> D_now; // current dye field
    std::vector<float> D_prev; // previous dye field 

    // frequency-space buffers
    fftwf_complex* u_freq;
    fftwf_complex* v_freq;

    // FFTW plans
    fftwf_plan u_plan_fwd, u_plan_bwd;
    fftwf_plan v_plan_fwd, v_plan_bwd;

    // counter for the number of simulation steps
    int count = 1; 

    void writeLog(std::ofstream& file, std::string logMessage) {
        file << logMessage << std::endl;
    }

    // Helper functions for solver member function (all logic from Jos Stam's algorithm):

    /* First, add external forces (if applicable) to the velocity fields.
    Then, update previous velocity fields to equal the new state*/
    void addForces() {
        for (int i = 0; i < n * n; i++) { 
            u_now[i] += dt * f_x[i];
            u_prev[i] = u_now[i];
            v_now[i] += dt * f_y[i];
            v_prev[i] = v_now[i];
            D_prev[i] = D_now[i];   
        }
    }

    /* Advect the velocity fields: semi-Lagrangian, stable even with large time step.
    1. for each voxel of the grid, calculate the new velocity u and v at time t
    2. trace voxel midppoint back to the previous time step t-dt (u_prev, v_prev)
    3. linearly interpolate the velocity at this point using neighbouring voxels
    4. assign the interpolated value to the voxel of step 1.
    */
    void advect() {
        int i, j, i0, j0, i1, j1;
        float x, y, x0, y0, s, t;
        for (x = 0.5 / n, i = 0; i < n; i++, x += 1.0 / n) {
            for (y = 0.5 / n, j = 0; j < n; j++, y += 1.0 / n) {

                // interpolate velocity at the voxel midpoint
                // x0 and y0 are the coordinates of the voxel in the previous time step
                x0 = n * (x - dt * u_prev[i + n * j]) - 0.5;
                y0 = n * (y - dt * v_prev[i + n * j]) - 0.5;
                i0 = floor(x0);
                s = x0 - i0;
                i0 = (n + (i0 % n)) % n;    // enforce periodic boundary conditions
                i1 = (i0 + 1) % n;

                j0 = floor(y0);
                t = y0 - j0;
                j0 = (n + (j0 % n)) % n;    
                j1 = (j0 + 1) % n;

                // linear interpolation of velocity fields
                u_now[i + n * j] = (1 - s) * ((1 - t) * u_prev[i0 + n * j0] + t * u_prev[i0 + n * j1]) +
                            s * ((1 - t) * u_prev[i1 + n * j0] + t * u_prev[i1 + n * j1]);

                v_now[i + n * j] = (1 - s) * ((1 - t) * v_prev[i0 + n * j0] + t * v_prev[i0 + n * j1]) +
                            s * ((1 - t) * v_prev[i1 + n * j0] + t * v_prev[i1 + n * j1]);

                // dye advection
                D_now[i + n * j] = (1 - s) * ((1 - t) * D_prev[i0 + n * j0] + t * D_prev[i0 + n * j1]) +
                             s * ((1 - t) * D_prev[i1 + n * j0] + t * D_prev[i1 + n * j1]);
            }
        }
        // for (i = 0; i < n; i++)
        // for (j = 0; j < n; j++)
        // {
        //     D_prev[i + n * j] = D_now[i + n * j];
        // }
    }

    void copyToFFT(float* spatial_data, float* fft) {
        // copy real data to FFTW buffer
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fft[i + n * j] = spatial_data[i + n * j];
            }
        }
    }

    enum class FFTDirection { Forward = 1, Inverse = -1 };

    // forward FFT: input: real array, output: complex array
    // inverse FFT: input: complex array, output: real array
    void performFFT(FFTDirection direction, float* real_data, fftwf_complex* complex_data, fftwf_plan plan) {
        if (direction == FFTDirection::Forward) {
            fftwf_execute_dft_r2c(plan, real_data, complex_data);
        } else {
            fftwf_execute_dft_c2r(plan, complex_data, real_data);
        }
    }

    /* Apply viscosity and project velocity field to conserve mass
    These are linear operations in the frequency domain.
    note: u_freq and v_freq are type fftwf_complex arrays
    1. apply viscous diffusion by filtering high frequencies
    2. project velocity field to conserve mass (divergence zero) by Helmholz decomposition
    */
    void applyViscosityAndProject() {
        float U_real, U_imag, V_real, V_imag;
        float x, r, f;
        for (int i = 0; i < (n/2+1); i++) {     // different loop than Stam's code due to different data structure
            x = i;                             // no *0.5 as we loop i++ and not i+=2
            for (int j = 0; j < n; j++) {
                int l = j;
                if (j > n / 2)      
                    l = j - n;

                r = x * x + l * l;
                if (r == 0)     // avoid division by zero
                    continue;

                // viscous diffusion
                f = expf(-r * dt * visc);

                // copy complex data to real and imaginary parts for calculations
                U_real = u_freq[i + (n/2+1) * j][0];
                U_imag = u_freq[i + (n/2+1) * j][1];
                V_real = v_freq[i + (n/2+1) * j][0];
                V_imag = v_freq[i + (n/2+1) * j][1];

                // apply mass conservation
                u_freq[i + (n/2+1) * j][0] = f * ((1 - x * x / r) * U_real - x * l / r * V_real);
                u_freq[i + (n/2+1) * j][1] = f * ((1 - x * x / r) * U_imag - x * l / r * V_imag);
                v_freq[i + (n/2+1) * j][0] = f * (-l * x / r * U_real + (1 - l * l / r) * V_real);
                v_freq[i + (n/2+1) * j][1] = f * (-l * x / r * U_imag + (1 - l * l / r) * V_imag);
            }
        }
    }

    // normalize result after inverse FFT
    void normalize() {
        float f;
        f = 1.0 / (n * n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                u_now[i + n * j] = f * u_fft[i + n * j];
                v_now[i + n * j] = f * v_fft[i + n * j];
                
            }
    }

 };

int main() {
    // simulation parameters
    const int n = 500;                       // grid size (n x n)
    const float viscosity = 0.001;           // viscosity
    const float dt = 0.01;                   // time step
    float time = 0.0f;
    int step = 0;
    const int maxSteps = 500;
    const int forceSteps = 10; // number of steps to apply external forces
    
    // force field parameters
    const float U0 = 5.0f;
    const float delta = 0.025f;
    const float A = 1.0f;
    const float sigma = 0.02f;
    const float k = 4.0f;

    // output file and log file parameters & initialisation
    int outputInterval = 10;
    int logInterval = 10;
    std::ofstream logFile;
    std::string logFilename = "log.txt";
    logFile.open(logFilename);
    logFile << "FluidSolver log" << std::endl;
    logFile << "Total simulation time: " << maxSteps * dt << " seconds" << std::endl;
    logFile << "Time step: " << dt << " seconds" << std::endl;
    logFile << "Grid size: " << n << " x " << n << std::endl;
    logFile << "Viscosity: " << viscosity << std::endl;
    std::ofstream outputFile;
    std::string filename = "output.csv";
    outputFile.open(filename);
    outputFile << "time;step;i;j;u;v;D" << std::endl;       
     

    //--- the main simulation loop ---
    FluidSolver solver(n, viscosity, dt, logFile);
    
    solver.setInitialDyeField(); 

    // first, run the simulation with external forces applied
    // note the mode of setting external forces (0 = set to zero, 1 = apply force from assignment)
    solver.setExternalForces(1, U0, delta, A, sigma, k);
    while (step < forceSteps) {
        solver.solve();
        
        if (step % outputInterval == 0) {
            solver.writeToFile(outputFile, time, step);
        }
        time += dt;
        step++;
    }
    // then, run the simulation without external forces
    solver.setExternalForces(0, U0, delta, A, sigma, k);
    while (step < maxSteps) {
        solver.solve();
        if (step % outputInterval == 0) {
            solver.writeToFile(outputFile, time, step);
        }
        time += dt;
        step++;
    }

    outputFile.close();
    logFile.close();
    std::cout << "Simulation completed. Output written to " << filename << std::endl;
    std::cout << "Log file written to " << logFilename << std::endl;
    return 0;
}
