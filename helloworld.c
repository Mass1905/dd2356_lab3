#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int rank;
    int nr_of_processes;
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nr_of_processes);
    
    printf("Hello world from rank %d from %d processes!\n", rank, nr_of_processes);
    MPI_Finalize();
}
