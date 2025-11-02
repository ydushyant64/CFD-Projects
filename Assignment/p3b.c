#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 41            // number of grid lines
#define mew 0.001       // dynamic viscosity
#define dx 1.0            // grid space
#define t_max 18        // maxmum time

void initialize(double x[], double u[])
{
    FILE *fp2 = fopen("Initial_values.txt","w");
    fprintf(fp2,"x \t\t\t u\n");
    //initial condition
    for(int i=0; i< N; i++)
    {
        x[i] = i*dx;
        u[i] = 0.5*(1 + tanh(250.0*(x[i] - 20.0)));
        fprintf(fp2,"%0.1lf \t %lf\n",x[i],u[i]);
    }
    fclose(fp2);
}

void boundary_condition(double u[])
{
     // initializing the (drichlet) boundary conditions
    u[0] = 0.0;            //u(0,t) = 0.0          
    u[N-1] = 1.0;          //u(40,t) = 1.0 
}

double mac_cormack(double u[], double u_next[], double dt)
{
    double u_bar[N];
    double dx2 = pow(dx,2.0);
    
    //Predictor term (forward differencing)
    for(int i=1; i<N -1; i++)
    {
        u_bar[i] = u[i] - ((dt/dx)*(0.5 - u[i])*(u[i+1] - u[i])) + ((dt*mew/dx2)*(u[i+1] - 2*u[i] + u[i-1]));
    }

    boundary_condition(u_bar);

    //Corrector term (backward differncing)
    for(int i = 1; i<N -1; i++)
    {
        u_next[i] = 0.5*(u[i] + u_bar[i] - (dt/dx)*((0.5 - u_bar[i])*(u_bar[i] - u_bar[i-1])) + (dt*mew/dx2)*(u_bar[i+1] - 2*u[i] + u[i-1]));
    }

    boundary_condition(u_next);
}

int main()
{
    double u[N], x[N], u_next[N];
    double dt = 0.5;
    double t=0.0;

    //Initializing the grid condiion
    initialize(x,u);

FILE *fp;
fp = fopen("delta_t_0.5.txt","w");

    while(t<t_max)
    {
        mac_cormack(u,u_next,dt);
        for(int i=0; i<N; i++)
        {
             u[i] = u_next[i];
        }

        t += dt;
    }
    printf("x \t\t\t u\n");
    fprintf(fp,"x \t\t\t u\n");
    for(int i=0; i<N; i++)
    {
        printf("%lf \t %lf\n",x[i],u[i]);
        fprintf(fp,"%lf \t %lf\n",x[i],u[i]);
    }

    fclose(fp);
    return 0;
}