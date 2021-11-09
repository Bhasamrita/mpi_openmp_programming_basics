/*The program adds the elements of an array in 3 ways-
  using block, cyclic, and block-cyclic distributions. 
  
  Different types of distribution of data is required because a particular
  code may contain dependency in different ways. For example, if a few consecutive
  lines are dependent on each other, then block distribution can be used. If every
  few blocks of data are dependent on each other, then we can use block-cyclic 
  distribution, and so on. */

#include <iostream>
#include <ctime>
#include <cmath>
#include <mpi.h>

using namespace std;

void share_local_sum (int &local_sum, int pid, int np);

int main(int argc, char** argv)
{
    srand((unsigned) time(0));
 
    MPI_Status status;
 
    MPI_Init(&argc, &argv);
 
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
 
    int arr_size;
    if (pid == 0)
    {
        //Initialise the array size (contains a maximum of 100 elements and minimum of 16)
        arr_size = rand() % 85 + 16;
     
        cout << "The array size is " << arr_size << endl;
    }
 
    MPI_Bcast(&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    int arr[arr_size];
    if (pid == 0)
    {
        for (int i = 0; i < arr_size; i++)
        {
            arr[i] = i;
            cout << arr[i] << " ";
        }
        cout << endl;
    }
 
    int size = ceil(arr_size / np);
    int block_array[size] = {0};
    
    //Block distribution of elements
    if (pid == 0)
    {
        cout << "The individual size each array gets is " << size << endl;
        for (int i = 0; i < np; i++)
        {
            int c = i * size;
            MPI_Send(&arr[c], size, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
     
    MPI_Recv(&block_array, size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    int local_sum = 0;
    for (int i = 0; i < size; i++)
    {
        cout << block_array[i] << " ";
        local_sum += block_array[i];
    }
    cout << endl;
    
    share_local_sum(local_sum, pid, np);
 
    if (pid == 0)
      cout << "For block distribution, the total sum of the array is - " << local_sum << endl;
 
    int cyclic_array[size] = {0};
 
    //Cyclic distribution of elements
    if (pid == 0)
    {
        for (int i = 0; i < np; i++)
        {
            int array[size] = {0};
            int k = 0;
            for (int j = i; j < arr_size; j += np)
                array[k++] = arr[j];
                
            MPI_Send(&array, size, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
    }
 
    //if (pid != np-1)
        MPI_Recv(&cyclic_array, size, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
 
    local_sum = 0;
    for (int i = 0; i < size; i++)
        local_sum += cyclic_array[i];
 
    share_local_sum(local_sum, pid, np);
 
    if (pid == 0)
      cout << "For cyclic distribution, the total sum of the array is - " << local_sum << endl;
 
    int sz, b;
    if (pid == 0)
    {
        b = ceil(size / np);
        sz = b * np;
    }
 
    MPI_Bcast(&sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    int bc_array[sz];

    //Block cyclic distribution of elements
    if (pid == 0)
    {
        for (int i = 0; i < np; i++)
        {
            int array[sz] = {0};
            int k = 0;
            for (int j = b * i; j < arr_size; j += ((b-1) * np))
            {
                int n = b;
                while(n-- && j < arr_size)
                    array[k++] = arr[j++];     
            }
        
            MPI_Send(&array, size, MPI_INT, i, 2, MPI_COMM_WORLD);
        }
    }
 
    MPI_Recv(&bc_array, size, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
 
    local_sum = 0;
    for (int i = 0; i < size; i++)
        local_sum += bc_array[i];
 
    share_local_sum(local_sum, pid, np);
 
    if (pid == 0)
      cout << "For block-cyclic distribution, the total sum of the array is - " << local_sum << endl;
 
    MPI_Finalize();
 
    return 0;
}

void share_local_sum(int &local_sum, int pid, int np)
{
    MPI_Status status;

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
}
