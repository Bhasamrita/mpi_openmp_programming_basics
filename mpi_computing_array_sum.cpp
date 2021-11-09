/*The following program computes the total size of an array in process 0 and then
  initialises the whole array in process 0. Then the array is divided into 
  similar sized blocks and each block is given to one of the processes. Each
  process then computes the sum of the elements in its block. The local_sum is
  then shared between processes in a tree fashion. The final sum is obtained in
  process 0, which then displays the output. */

#include <iostream>
#include <ctime>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    srand((unsigned) time(0));    

    int size;
 
    MPI_Status status;
 
    MPI_Init(&argc, &argv);
 
    //Initialising process id
    int pid;
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
 
    //Accepting the total number of processes to run
    int np;
    MPI_Comm_size (MPI_COMM_WORLD, &np);
 
    if (pid == 0)
    {
        //The size of array and individual elements are so chosen such that
        //the sum is contained in an integer.
     
        //Enter the number of elements in the array
        int arr_size = rand() % 500 + 1;
        cout << "Array size is " << arr_size << endl;
     
        size = ceil(arr_size / np);
        //cout << "Size = " << size << endl;
     
        //Broadcasting the size of array that each process will contain at most
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
     
        int array[arr_size];
        //Initialising the array elements
        for (int i = 0; i < arr_size; i++)
        {
            //if (i % 20 == 0)
              //cout << endl;
         
            array[i] = rand() % 100 + 1;
            //cout << array[i] << " ";
        }
        //cout << endl;
     
        //A chunk of array is sent to every process
        for (int i = 0; i < np; i++)
        {
            int c = i * size;
            MPI_Send(&array[c], size, MPI_INT, i, 1, MPI_COMM_WORLD);
            //cout << "I have sent!" << endl;
        }
    }
    else
    {
        //Every process receives the maximum size of array it will contain
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //cout << "Process " << pid << " received size = " << size << endl;
    }
  
    //Every process receives the array elements
    int arr[size] = {0};
    MPI_Recv(&arr, size, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    cout << "Process " << pid << " received array." << endl;
 
    //Every process calculates the local_sum of its elements
    int local_sum = 0;
    for (int i = 0; i < size; i++)
        local_sum += arr[i];
    
    //The local_sum is shared between processes in a tree fashion
    int j = 1;
    while ((int)pow(2, j) <= np)
    {
        //cout << "pid = " << pid << " j = " << j << endl;
        if (pid % (int)pow(2,j) == 0)
        {
            /*The receiving process initialises a buffer to receive the data
              and then adds it to its local sum variable. */
            int temp;
            MPI_Recv (&temp, 1, MPI_INT, pid+j, 0, MPI_COMM_WORLD, &status);
            local_sum += temp;
            //cout << "pid = " << pid << " local_sum = " << local_sum << endl;
        }
        else if (pid % (int)pow(2,j-1) == 0)
        {
            MPI_Send (&local_sum, 1, MPI_INT, pid-j, 0, MPI_COMM_WORLD);
            //cout << pid << " have sent!" << endl;
            break;
        }
        j++;
    }
 
    //The last process contains the total sum and it displays the same
    if (pid == 0)
        cout << "The total sum is: " << local_sum << endl;
 
    MPI_Finalize();
 
    return 0;
}
