#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 21
#define C0 0.4
#define r 0.1
#define max_iter 1000
#define epsilon 1e-6

double initialization(double u[],double u_new[]);
double analytical(double u[],double dx);
double FTCS(double u[], double u_new[]);
double MacCormack(double u[], double u_new[]);


//FTCS Scheme
double FTCS(double u[], double u_new[N])
{
    double error[max_iter];
    int iter;
    for(iter = 0; iter<max_iter; iter++)
    {
    for(int i =1; i<N-1; i++)
    {
        u_new[i] = u[i] - 0.5*C0*(u[i+1] - u[i-1]) + r*(u[i+1] -2*u[i] + u[i-1]);   
        
        for(int j = 1;j <N-1; j++)
            error[iter] += pow(u_new[j] - u[j],2);
        error[iter] = pow(error[iter],0.5);
        for(int j = 1; j < N-1; j++)
            u[j] = u_new[j];
        if(error[iter] < epsilon)
            break;
    }
    } 
    return iter;
}

//Mac-Cormack method
double MacCormack(double u[], double u_new[])
{   
   double u_bar[N], error[max_iter];
   int iter;
    u_bar[0] = u[0];
    u_bar[N-1] = u[N-1];

    for(iter = 0; iter < max_iter; iter++)
    {
        //Predictor term (forward differncing)
        for(int i = 1; i < N-1; i++)
        {
            u_bar[i] = u[i] - C0*(u[i+1] - u[i]) + r*(u[i+1] -2*u[i] + u[i-1]);
        }


        //Corrector term (Backward differencing)
        for(int i = 1; i < N-1; i++)
        {
            u_new[i] = 0.5*(u[i] + u_bar[i] - C0*(u_bar[i] - u_bar[i-1]) + r*(u_bar[i+1] - 2*u_bar[i] + u_bar[i-1]));
        }


        for(int j = 1;j <N-1; j++)
            error[iter] += pow(u_new[j] - u[j],2);
        error[iter] = pow(error[iter],0.5);
        for(int j = 1; j < N-1; j++)
            u[j] = u_new[j];
        if(error[iter] < epsilon)
            break;
    }
    return iter;
}

double analytical(double u[],double dx)
{
    double x[N];
    for (int i=0; i<N; i++)
    {
        x[i] = i*dx;
        u[i] = 100.0*((1 - exp((C0/(r*dx))*(x[i] - 1)))/(1 - exp(-C0/(r*dx))));
    }
}


// Initializing the boundary conditions
double initialization(double u[],double u_new[])
{
    for(int i=0; i<N; i++)
    {
        if(i ==0)
        {
          u[i] = 100.0;
          u_new[i] = 100.0;
        }
        else
        {
           u[i] = 0.0;
           u_new[i] = 0.0;
        }
    }
    return 0;
}

int main()
{   
    double U_ftcs[N],U_mc[N],U_an[21],U_new[N],x[N];
    double dx = 1.0/20.0;
    
    initialization(U_ftcs,U_new);
    int iter = FTCS(U_ftcs,U_new);    
    printf("FTCS scheme converges in %d iterations\n",iter);

    initialization(U_mc,U_new);
    iter = MacCormack(U_mc,U_new);
    printf("MacCormack scheme converges in %d iterations\n",iter);

    analytical(U_an,dx);
    
    FILE *fp = fopen("Output2.txt","w");
    fprintf(fp,"x \t\t U_FTCS \t U_MacCormack \t U_Analytical\n");
    for(int i = 0; i < N; i++)
    {
        x[i] = i*dx;
        fprintf(fp,"%0.4lf \t %0.4lf \t %0.4lf \t %0.4lf\n",x[i],U_ftcs[i],U_mc[i],U_an[i]);
    }
    fclose(fp);
    return 0;
}

