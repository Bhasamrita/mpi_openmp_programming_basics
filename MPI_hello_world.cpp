#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    // pid is process id and np is the total number of processes
    int pid, np;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
    MPI_Comm_size (MPI_COMM_WORLD, &np);

    cout << "Hello World!";
    cout << "I am process " << pid << " out of " << np << " processes." << endl;

    MPI_Finalize();
}
