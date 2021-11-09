//This program tries to duplicate the working of MPI_Barrier

#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    MPI_Status status;
 
    MPI_Init (&argc, &argv);
 
    int pid;
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
 
    int np;
    MPI_Comm_size (MPI_COMM_WORLD, &np);
 
    cout << "I am process " << pid << " and I am about to enter the barrier." << endl;
 
    //Implementing iittp_barrier
    char should_I_get_out;
    int I_am_in;
    if (pid == 0)
    {
        int count = 1;
        while (count < np)
        {
            MPI_Recv(&I_am_in, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            count ++;
            //cout << count << " ";
        }
        cout << endl;

        should_I_get_out = 'y';
        MPI_Bcast(&should_I_get_out, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    else
    {
        I_am_in = 1;
        MPI_Send(&I_am_in, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
         
        MPI_Bcast(&should_I_get_out, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    }

    cout << "I am process " << pid << " and I am out of the barrier!" << endl;
 
    MPI_Finalize();
 
    return 0;
}
