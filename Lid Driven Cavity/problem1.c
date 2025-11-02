#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Constants
#define Nx 32
#define Ny Nx
#define L 1.0
#define Wall_Velocity 4.0
#define rho 1.0
#define mu 0.01
#define dt 0.001
#define maxIt 50000
#define maxe 1e-7

// Helper function to initialize arrays
void initialize(double arr[Nx][Ny], double value) 
{
    for (int i = 0; i < Nx; i++) 
    {
        for (int j = 0; j < Ny; j++) 
        {
            arr[i][j] = value;
        }
    }
}

// Main solver
int main() 
{
    // Variables
    double Vo[Nx][Ny], St[Nx][Ny], Vop[Nx][Ny], u[Nx][Ny], v[Nx][Ny];
    double h = (double)L / (Nx - 1);
    double error = 0.0;

    // Initialize arrays
    initialize(Vo, 0.0);
    initialize(St, 0.0);
    initialize(Vop, 0.0);
    initialize(u, 0.0);
    initialize(v, 0.0);

    // Iterative solver (Gauss-Seidel-like)
    for (int iter = 0; iter < maxIt; iter++) 
    {
        // Apply boundary conditions
        for (int i = 0; i < Nx; i++) 
        {
            Vo[i][Ny - 1] = -2.0 * St[i][Ny - 2] / (h * h) - Wall_Velocity * 2.0 / h; // Top
            Vo[i][0] = -2.0 * St[i][1] / (h * h);                                   // Bottom
        }
        for (int j = 0; j < Ny; j++) 
        {
            Vo[0][j] = -2.0 * St[1][j] / (h * h);                                   // Left
            Vo[Nx - 1][j] = -2.0 * St[Nx - 2][j] / (h * h);                         // Right
        }

        // Update vorticity using transport equation
        for (int i = 1; i < Nx - 1; i++) 
        {
            for (int j = 1; j < Ny - 1; j++) 
            {
                Vop[i][j] = Vo[i][j];
                Vo[i][j] = Vop[i][j] +(-((St[i][j + 1] - St[i][j - 1]) / (2 * h)) *((Vop[i + 1][j] - Vop[i - 1][j]) / (2 * h)) +((St[i + 1][j] - St[i - 1][j]) / (2 * h)) *((Vop[i][j + 1] - Vop[i][j - 1]) / (2 * h)) +mu / rho *((Vop[i + 1][j] + Vop[i - 1][j] - 4 * Vop[i][j] +Vop[i][j + 1] + Vop[i][j - 1]) /(h * h))) *dt;
            }
        }

        // Update stream function using elliptical equation
        for (int i = 1; i < Nx - 1; i++) 
        {
            for (int j = 1; j < Ny - 1; j++) 
            {
                St[i][j] = (Vo[i][j] * h * h + St[i + 1][j] + St[i - 1][j] +St[i][j + 1] + St[i][j - 1]) /4.0;
            }
        }

        // Check convergence
        if (iter > 10)
        {
            error = 0.0;
            for (int i = 1; i < Nx - 1; i++) 
            {
                for (int j = 1; j < Ny - 1; j++) 
                {
                    double diff = fabs(Vo[i][j] - Vop[i][j]);
                    if (diff > error) 
                    {
                        error = diff;
                    }
                }
            }
            if (error < maxe) 
            {
                printf("For Re = 400, Vorticity Converged at iteration %d with error %e\n", iter, error);
                break;
            }
        }
    }

    // Calculate velocities
    for (int i = 1; i < Nx - 1; i++) 
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            u[i][j] = (St[i][j + 1] - St[i][j - 1]) / (2 * h);
            v[i][j] = -(St[i + 1][j] - St[i - 1][j]) / (2 * h);
        }
    }

    // Output results to file (for visualization)
    FILE *output = fopen("results.txt", "w");
    fprintf(output, "x \t\t\t y \t\t\t u \t\t\t v \t\t\t Vorticity\tStream_function\n");
    for (int i = 0; i < Nx; i++) 
    {
        for (int j = 0; j < Ny; j++) 
        {
            fprintf(output, "%f\t%f\t%f\t%f\t%f\t%f\n", i * h, j * h, u[i][j], v[i][j], Vo[i][j], St[i][j]);
        }
    }
    fclose(output);
    return 0;
}
