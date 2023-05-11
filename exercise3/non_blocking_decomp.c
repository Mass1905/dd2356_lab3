#include <mpi.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define  MOD(a, b)  ((((a)%(b))+(b))%(b))




int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 900; // make sure nxc is divisible by size
    double L = 2*3.141592653589793; // Length of the domain
    

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nxc % size){
      printf(" Error, nxc must be devisible by size \n");
      return (1);
    }
    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++){
      f[i] = sin(L_loc*rank + (i-1) * dx); 
    }
    /*
    if(rank == size -1)
      f[nxn_loc-2]=0;
  */
    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
    MPI_Status statuses [2];
    MPI_Request requests_recv[2];
    MPI_Request requests_send[2];

    //send msgs 
    MPI_Isend(f+2, 1, MPI_DOUBLE, MOD(rank - 1, size), 0, MPI_COMM_WORLD, requests_send);
    MPI_Isend(f + nxn_loc - 3, 1, MPI_DOUBLE, MOD(rank + 1, size), 1, MPI_COMM_WORLD, requests_send + 1);

    //receive msgs
    MPI_Irecv(f + nxn_loc - 1, 1, MPI_DOUBLE, MOD(rank + 1, size), 0, MPI_COMM_WORLD, requests_recv);
    MPI_Irecv(f, 1, MPI_DOUBLE,  MOD(rank - 1, size), 1, MPI_COMM_WORLD, requests_recv + 1);

    //wait until ghost cells values are received
    MPI_Waitall(2, requests_recv, statuses );

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);


    sleep(rank + 3);
    for (i=1; i< (nxn_loc - 1 ) ; i++){
        double cc = cos(L_loc*rank +  (i-1) * dx);
        double diff = (cc - dfdx[i] >= 0)? (cc - dfdx[i]) : ( dfdx[i] - cc);
        printf("sin'(%f) == %f \n", L_loc*rank + (i-1) * dx, dfdx[i]);

        if( diff > (1 / (double)nxc ))
          printf("difference detected at x= %f \n", L_loc*rank +  (i-1) * dx);
  
       }
    MPI_Finalize();
}






