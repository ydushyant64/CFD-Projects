#include <stdio.h>
#include <math.h>
#define rho 1.2 //for air (kg/m3).
#define P_atm 0 // taking pressure at A-B as atm pressure.
#define U0 1 //free stream velocity. 
#define l_x 10
#define l_y 10
#define dx 0.5
#define dy 0.5



void point_gs(int nx, int ny, double psi[nx][ny], double epsilon);

int main() 
{
    int nx = (int)(l_x / dx) + 1;
    int ny = (int)(l_y / dy) + 1;
    double psi[nx][ny];
    double epsilon = 1e-7;

    // Initialize the stream function array to zero
    for (int i = 0; i < nx; i++) 
    {
        for (int j = 0; j < ny; j++) 
        {
            psi[i][j] = 0;
        }
    }

    // Set boundary conditions
    for (int i = 0; i < nx; i++) 
    {
        psi[i][0] = 0;       // Bottom boundary (y = 0, ψ = 0)
        psi[i][ny - 1] = 10; // Top boundary (y = H, ψ = 10)
    }

    // Inlet (A-B) boundary condition with linear stream function increase
    for (int j = 0; j < ny; j++) 
    {
        psi[0][j] = U0 * j * dy;  // Left boundary (A-B)
    }

    point_gs(nx, ny, psi, epsilon);

    return 0;
}

void point_gs(int nx, int ny, double psi[nx][ny], double epsilon) 
{
    int iter = 0;
    double psi_old, psi_new;
    double sum;
    int x_center = nx - 1;  // Bottom-right corner
    int y_center = 0;
    double radius = 2.0;
    double v[ny];
    double p;
    double cp;

    FILE *fp1, *fp2, *fp3, *fp4;

    // Open files to store results
    fp1 = fopen("point_gs_convergence.txt", "w");
    fp2 = fopen("point_gs_streamfunction.txt", "w");
    fp3 = fopen("velocity_at_C-D.txt", "w");
    fp4 = fopen("coefficient_of_pressure_and_comparison.txt", "w");
     // Iterative Gauss-Seidel loop
    do {
        sum = 0;
        for (int i = 1; i < nx - 1; i++) 
        {
            for (int j = 1; j < ny - 1; j++) 
            {
                // Mask region inside a circle of radius 2 centered at (x_center, y_center)
                double x_dist = (i - x_center) * dx;
                double y_dist = (j - y_center) * dy;

                if (x_dist * x_dist + y_dist * y_dist >= radius * radius) 
                {
                psi_old = psi[i][j];
                psi_new = 0.25 * (psi[i - 1][j] + psi[i + 1][j] + psi[i][j - 1] + psi[i][j + 1]);
                double diff = psi_old - psi_new;
                psi[i][j] = psi_new;
                sum += diff * diff;  // Accumulate squared error
                psi[nx-1][j] = psi[nx-2][j];
                }
            }
        }

        sum = sqrt(sum);  // Root of the sum of squared errors
        fprintf(fp1, "%d\t%.8lf\n", iter, sum);  // Write convergence data
        iter++;
    } 
    while (sum > epsilon);
    fclose(fp1);

    // Write final stream function distribution
    printf("Converged after %d iterations with maximum error = %e\n", iter, sum);
    printf("Stream function distribution:\n");
    for (int i = 0; i < nx; i++) {  // Top to bottom for visualization
        for (int j = 0; j < ny; j++) {
            printf("%6.2f ", psi[j][i]);
            fprintf(fp2, "%6.2f ", psi[j][i]);
        }
        printf("\n");
        fprintf(fp2, "\n");
    }
    //To calculate the velocity at C-D
    for (int i = 0; i < ny-1; i++)
    {
        v[0]=0;
        v[i+1] = (psi[nx-1][i+1] - psi[nx-1][i])/dy;
        fprintf(fp3, "velocity[%d]: %6.2f\n", i, v[i]);
        if (dy*i == radius)
        {
            p = P_atm + 0.5 * rho * (pow(U0,2) - v[i]*v[i]);
            cp = (p - P_atm) / (0.5*rho*pow(U0,2));
            fprintf(fp4, "pressure at point D: %6.2lf kpa\n", p);
            fprintf(fp4, "coefficient of pressure at point D: %6.2lf\n", cp);
        }
        
    }
    
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
}
