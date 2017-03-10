#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define PI acos(-1.0)

// Parameters

int N = 64;
double delta_t  = 5E-3;
double betta = 1.0;

// Define internal functions

double* init(double *array);
double omega(int k);

// Main program

int main(int argc, char *argv[]){

    clock_t tic = clock();
  	int T = 1.0*pow(N,2.2);
    int i, j, k, t, n, cont = 0, cont2 = 0;
    int P = 1000, K = 3, V = 200;
    double Q, Qp;
    int iterations = (int) T/delta_t;

    // Setting the number of threads
    int threads = atoi(argv[1]);
    omp_set_num_threads(threads);
    printf("we are running with %d threads\n", threads);

  	// Create the grids to store N points in the solid

  	double* x = malloc(N*sizeof(double));
  	double* x_new = malloc(N*sizeof(double));
  	double* v = malloc(N*sizeof(double));
  	double* v_new = malloc(N*sizeof(double));
  	double* F_grid = malloc(N*sizeof(double));
  	double* F_grid_new = malloc(N*sizeof(double));

  	// Create the total grid to store data

    double **energy = (double **) malloc(K * sizeof(double *));
    for(i=0; i<K; i++){
		    energy[i] = (double *) malloc(P * sizeof(double *));
	  }
    double **DATA = (double **) malloc(V * sizeof(double *));
    for(i=0; i<V; i++){
        DATA[i] = (double *) malloc(N * sizeof(double *));
    }

  	// Initial conditions for the solid
    // Position
	  x = init(x);
  	x_new = init(x_new);
    // Velocity & Acceleration
    for(n=0;n<N;n++){
      v[n] = 0.0;
      v_new[n] = 0.0;
      F_grid[n] = 0.0;
      F_grid_new[n] =0.0;
    }

    // ---------------- Leapfrog method ---------------------- //
	  // We use the Leapfrog method expressed in x, v & a quantities with integer steps
    for(t = 0; t < iterations; t++){

      #pragma omp parallel for private(n), shared(F_grid,F_grid_new,x,x_new,v,v_new,delta_t,N)
    	for(n = 1; n < N-1; n++){
    		F_grid[n] = (x[n+1] - 2.0 * x[n] + x[n-1]) + betta * (pow((x[n+1] - x[n]),2.0) - pow((x[n] - x[n-1]),2.0) );
    	  x_new[n] = x[n] + ( v[n] * delta_t ) + (0.5 * F_grid[n] * pow(delta_t,2.0) );
        F_grid_new[n] = (x_new[n+1] - 2.0 * x_new[n] + x_new[n-1]) + betta * (pow((x_new[n+1] - x_new[n]),2.0) - pow((x_new[n] - x_new[n-1]),2.0));
        v_new[n] = v[n] + 0.5 * (F_grid[n] + F_grid_new[n]) * delta_t;
      }

      // Update
      x = x_new;
      v = v_new;

      // -------------------- Store Data ----------------------- //
      // Print the velocities and positions
      if(t%(iterations/P) == 0){
        // Print iterations
//        printf("iter: %d  from %d , cont :%d\n",t,iterations, cont);

        // Calculate the general position
        Q = 0.0; Qp = 0.0;
//        #pragma omp parallel for private(n)
        for ( k = 1; k <= K; k++) {
          for ( n = 0; n < N; n++) {
            //            printf("%d %d\n",k, n);
            Q = Q + x[n]*sin((PI*(double)k*n)/(double)(N+1));
            Qp = Qp + v[n]*sin((PI*(double)k*n)/(double)(N+1));
          }
          Q = Q * sqrt(2.0/(double)(N+1));
          Qp = Qp * sqrt(2.0/(double)(N+1));

          energy[k-1][cont] = 0.5 * (pow(Q,2.0)+pow(omega(k)*Qp,2.0));
//          printf("%f %f\n", energy[k-1][cont], Qp );

        }
        cont++;
      }
      // Position
      if(t%(iterations/V) == 0 && cont2<V){
//        printf("Saving data position, cont2 : %d\n", cont2);
        for (int i = 0; i < N; i++) {
          DATA[cont2][i] = x[i];
        }
        cont2++;
      }
      // -------------------- Store Data ----------------------- //
    }
    // ---------------- Leapfrog method ---------------------- //

    // ------------------- Print results --------------------- //
    // Energy
    FILE * file1 = fopen("energy.dat","w");
    for(i=0; i<K ; i++){
      for(j = 0; j<P;j++){
        fprintf(file1,"%f \n", energy[i][j]);
      }
    }
    fclose(file1);

    // Position
    FILE * file2 = fopen("position.dat","w");
    for(i=0; i<V ; i++){
      for(j = 0; j<N;j++){
        fprintf(file2,"%f \n", DATA[i][j]);
      }
    }
    fclose(file2);
    // ------------------- Print results --------------------- //
    clock_t toc = clock();
    printf("%f\n", (double)(toc - tic) / CLOCKS_PER_SEC);
  	return(0);
}

// Auxiliar functions

double *init(double *array){
  	array[0] = 0.0;
  	array[N-1] = 0.0;
	  for(int n=1; n<N-1; n++){
		    array[n] = sin(PI*(double)n/(double)(N-1));
	  }
	  return array;
}
double omega(int k){
  return 2*sin(PI*k/((double)((2*N)+2)));
}
