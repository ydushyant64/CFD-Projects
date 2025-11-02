#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GRID_X 100      // Number of grid points in x
#define GRID_Y 50       // Number of grid points in y
#define LENGTH_X 2.0    // Domain length in x
#define LENGTH_Y 1.0    // Domain length in y
#define REYNOLDS 100.0  // Reynolds number
#define TIME_STEP 0.01  // Time step size
#define MAX_ITERATIONS 5000 // Maximum iterations
#define CONV_TOLERANCE 1e-6 // Convergence tolerance

void initialize_fields(double velocity_u[GRID_Y + 2][GRID_X + 1], double velocity_v[GRID_Y + 1][GRID_X + 2], double pressure[GRID_Y + 1][GRID_X + 1]) 
{
    for (int i = 0; i < GRID_Y + 2; i++) 
    {
        for (int j = 0; j < GRID_X + 1; j++) 
        {
            velocity_u[i][j] = 0.0;
        }
    }
    for (int i = 0; i < GRID_Y + 1; i++) 
    {
        for (int j = 0; j < GRID_X + 2; j++) 
        {
            velocity_v[i][j] = 0.0;
        }
    }
    for (int i = 0; i < GRID_Y + 1; i++) 
    {
        for (int j = 0; j < GRID_X + 1; j++) 
        {
            pressure[i][j] = 0.0;
        }
    }
    // Boundary conditions for inlet
    for (int i = 0; i < GRID_Y + 2; i++) 
    {
        velocity_u[i][0] = 1.0;  // Inlet velocity
    }
}

void solve_pressure_poisson(double pressure[GRID_Y + 1][GRID_X + 1], double rhs[GRID_Y + 1][GRID_X + 1], double dx, double dy) 
{
    for (int k = 0; k < 100; k++) 
    {
        for (int i = 1; i < GRID_Y; i++) 
        {
            for (int j = 1; j < GRID_X; j++) 
            {
                pressure[i][j] = 0.25 * (pressure[i + 1][j] + pressure[i - 1][j] +pressure[i][j + 1] + pressure[i][j - 1] -rhs[i][j] * dx * dy);
            }
        }
    }
}

double compute_max_divergence(double rhs[GRID_Y + 1][GRID_X + 1]) 
{
    double max_div = 0.0;
    for (int i = 1; i < GRID_Y; i++) 
    {
        for (int j = 1; j < GRID_X; j++) 
        {
            if (fabs(rhs[i][j]) > max_div) 
            {
                max_div = fabs(rhs[i][j]);
            }
        }
    }
    return max_div;
}

int main() 
{
    double dx = LENGTH_X / GRID_X, dy = LENGTH_Y / GRID_Y;
    double velocity_u[GRID_Y + 2][GRID_X + 1], velocity_v[GRID_Y + 1][GRID_X + 2], pressure[GRID_Y + 1][GRID_X + 1];
    double velocity_u_star[GRID_Y + 2][GRID_X + 1], velocity_v_star[GRID_Y + 1][GRID_X + 2], rhs[GRID_Y + 1][GRID_X + 1];

    initialize_fields(velocity_u, velocity_v, pressure);

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) 
    {
        // Predictor step
        for (int i = 1; i < GRID_Y + 1; i++) 
        {
            for (int j = 1; j < GRID_X; j++) 
            {
                velocity_u_star[i][j] = velocity_u[i][j] -TIME_STEP * ((velocity_u[i][j] * (velocity_u[i][j + 1] - velocity_u[i][j - 1]) / (2.0 * dx)) +(velocity_v[i][j] * (velocity_u[i + 1][j] - velocity_u[i - 1][j]) / (2.0 * dy)) -(1.0 / REYNOLDS) * ((velocity_u[i][j + 1] - 2.0 * velocity_u[i][j] + velocity_u[i][j - 1]) / (dx * dx) +(velocity_u[i + 1][j] - 2.0 * velocity_u[i][j] + velocity_u[i - 1][j]) / (dy * dy)));
            }
        }

        for (int i = 1; i < GRID_Y; i++)
        {
            for (int j = 1; j < GRID_X + 1; j++) 
            {
                velocity_v_star[i][j] = velocity_v[i][j] -TIME_STEP * ((velocity_u[i][j] * (velocity_v[i][j + 1] - velocity_v[i][j - 1]) / (2.0 * dx)) +(velocity_v[i][j] * (velocity_v[i + 1][j] - velocity_v[i - 1][j]) / (2.0 * dy)) -(1.0 / REYNOLDS) * ((velocity_v[i][j + 1] - 2.0 * velocity_v[i][j] + velocity_v[i][j - 1]) / (dx * dx) +(velocity_v[i + 1][j] - 2.0 * velocity_v[i][j] + velocity_v[i - 1][j]) / (dy * dy)));
            }
        }

        // Compute RHS for Poisson equation
        for (int i = 1; i < GRID_Y; i++) 
        {
            for (int j = 1; j < GRID_X; j++) 
            {
                rhs[i][j] = ((velocity_u_star[i][j + 1] - velocity_u_star[i][j]) / dx +(velocity_v_star[i + 1][j] - velocity_v_star[i][j]) / dy) / TIME_STEP;
            }
        }

        // Solve Poisson equation
        solve_pressure_poisson(pressure, rhs, dx, dy);

        // Corrector step
        for (int i = 1; i < GRID_Y + 1; i++) 
        {
            for (int j = 1; j < GRID_X; j++) 
            {
                velocity_u[i][j] = velocity_u_star[i][j] - TIME_STEP * (pressure[i][j] - pressure[i][j - 1]) / dx;
            }
        }

        for (int i = 1; i < GRID_Y; i++) 
        {
            for (int j = 1; j < GRID_X + 1; j++) 
            {
                velocity_v[i][j] = velocity_v_star[i][j] - TIME_STEP * (pressure[i][j] - pressure[i - 1][j]) / dy;
            }
        }

        // Check for convergence
        double divergence = compute_max_divergence(rhs);
        if (divergence < CONV_TOLERANCE) 
        {
            printf("Converged in %d iterations\n", iteration);
            break;
        }
    }

    // Output velocity field
    FILE *output_file = fopen("velocity.dat", "w");
    for (int i = 0; i < GRID_Y + 2; i++) 
    {
        for (int j = 0; j < GRID_X + 1; j++) 
        {
            fprintf(output_file, "%f ", velocity_u[i][j]);
        }
        fprintf(output_file, "\n");
    }
    fclose(output_file);

    return 0;
}
