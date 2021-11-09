/*The program computes the length of the array in process 0. The array cannot be
  stored in a single process. The array length is shared to all processes. Each
  process takes in (array_len / np) number of elements and computes the local sum.
  The local_sum for all processes are then shared to finally arrive at total_sum. */

#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
    //Initialising rand() function to system time
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
        /*The total size of array is being initialised in process 0. For 
          simplicity, the array size is assumed to be between 11 and 100. 
          It also makes sure that the sum doesn't overflow in an int variable. */
        arr_size = rand() % 99 + 11;
        cout << "Array size is " << arr_size << endl;
        
        //The array size is broadcasted
        MPI_Bcast (&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //cout << "Max size of array for every process = " << c << endl;
    }
    else
    {
        //All other process receives the value of arr_size
        MPI_Bcast (&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
 
    vector<int> vec;
    int local_sum = 0, c = ceil(arr_size / np);
    if (pid != np - 1)
    {
        //The array elements are initialised with the index of the complete array.
        for (int i = pid * c; i < (pid + 1) * c; i++)
            vec.push_back(i);
     
        /*Every process then adds up its array elements and stores the sum in 
          a local sum variable. */
        for (auto x : vec)
          local_sum += x;
        //cout << "pid = " << pid << " sum = " << local_sum << endl;
    }
    else
    {
        //The last process may contain less than c elements.
        for (int i = pid * c; i < arr_size; i++)
            vec.push_back(i);
     
        for (auto x : vec)
          local_sum += x;
        //cout << "pid = " << pid << " sum = " << local_sum << endl;
    }

    /*If we have 8 processes, then the below code ensures the following -
      - In the first iteration, all even processes will receive, 
        and odd processes will send and exit the loop. 
      - In the second iteration, all 4modulo processes will receive,
        and the remaining will send and exit.
      - In the third iteration, all 8modulo processes will receive, 
        and the remaining will send and exit.
      This loop will run for log_base_2(np) times. */
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
 
    /*In the end, the local sum variable of process 0 will contain the total
      sum of the array elements. */
    if (pid == 0)
      cout << "The total sum is: " << local_sum << endl;

    //Cleaning up all MPI states before exiting
    MPI_Finalize();
 
    return 0;
}
