#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int rank;
    int size;

    MPI_Init(NULL, NULL);
    double startTime, stopTime, elapsedTime;
    startTime = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int localCount = 0;
    int totalCount = 0;
    double x, y, z, pi;
    
    srand(SEED * rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < (NUM_ITER / size); iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            localCount++;
        }
    }

    MPI_Reduce(&localCount, &totalCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    stopTime = MPI_Wtime();
    elapsedTime = stopTime - startTime;

    if (rank == 0)
    {
        // Estimate Pi and display the result
        pi = ((double)totalCount / (double)(NUM_ITER)) * 4.0;

        printf("The result is %f\n", pi);
        printf("The time is %f\n", elapsedTime);
    }
    
    MPI_Finalize();
    return 0;
}