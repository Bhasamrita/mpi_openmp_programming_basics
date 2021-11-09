/*This process initialises an array size, which will be common to all the processes,
  in process 0 and then broadcasts it. Then, each process initialises that many
  number of elements. Then, process k sends its data to process k+1. The last 
  process sends its data to process 0. The array elements contained in all the
  processes are then displayed. */

#include <iostream>
#include <ctime>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    srand((unsigned) time(0));
 
    MPI_Status status;

    //Initialising the MPI environment
    MPI_Init (&argc, &argv);
 
    //Rank or process id
    int pid;
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
 
    //Size or number of processes
    int np;
    MPI_Comm_size (MPI_COMM_WORLD, &np);
 
    int arr_size;  
    if (pid == 0)
    {
        /*The total size of array is being initialised in process 0 and is
          broadcasted. I was having trouble taking input from the user so 
          I have used rand() function here. For simplicity, the array size 
          is assumed to be between 11 and 100. */
        arr_size = rand() % 90 + 11;
     
        MPI_Bcast (&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        cout << "Array size = " << arr_size << endl;
    }
    else
    {
        //All other process receives the value of arr_size
        MPI_Bcast (&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
 
    int local_arr[arr_size];

    /*Every process creates a local array of size arr_size and generates 
      that many random values. For simplicity again, the values are 
      considered to be between 1 and 10. */
      for (int i = 0; i < arr_size; i++)
          local_arr[i] = rand() % 10 + 1;

      //Process k sends its data to process k+1
      MPI_Send (&local_arr, arr_size, MPI_INT, (pid+1)%np, 1, MPI_COMM_WORLD);
      //cout << pid << " has sent data!" << endl;
 
      int temp[arr_size];
      //Process k receives its data from process k-1
      MPI_Recv (&temp, arr_size, MPI_INT, (pid+np-1)%np, 1, MPI_COMM_WORLD, &status);
      //cout << pid << " has received data!" << endl;
 
      //Displaying the arrays initialised by every process
      for (int i = 0; i < arr_size; i++)
        cout << local_arr[i] << " ";
      cout << endl;
 
      //Displaying the arrays received by every process
      for (int i = 0; i < arr_size; i++)
        cout << temp[i] << " ";
      cout << endl;
 
      //Cleaning up all MPI states before exiting
      MPI_Finalize();
 
      return 0;  
}
