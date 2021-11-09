/*This program initialises the vector size in a process and then sends the size
  to all other processes. The processes then calculate their share of the vector
  size and initialises the vector elements of the two vectors. The partial dot
  product is found out by all the processes and the sum is shared among the
  processes in a tree fashion. The final dot product is contained within process 0
  and it displays the output. */

#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
    srand((unsigned) time(0));
 
    MPI_Status status;
 
    MPI_Init(&argc, &argv);
 
    //Initialising process id
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 
    //Initialising total number of processes
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
 
    int vec_size;
    if (pid == 0)
    {
        //The vector size is initialised to a random number between 11 and 500
        vec_size = rand() % 490 + 11;
        cout << "The size of the vectors is " << vec_size << endl;
     
        //The vector size is broadcasted to all the processes
        MPI_Bcast(&vec_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
        //The vector size is received by all the processes
        MPI_Bcast(&vec_size, 1, MPI_INT, 0, MPI_COMM_WORLD);  
    }
 
    vector<int> vec1, vec2;
    int n = ceil(vec_size / np), dot_prod = 0;
 
    if (pid != np-1)
    {
        //The elements of the vector are kept small to avoid overflow of the dot product
        for (int i = 0; i < n; i++)
        {
            int val = rand() % 10;
            vec1.push_back(val);
         
            val = rand() % 10;
            vec2.push_back(val);
         
            dot_prod += (vec1[i] * vec2[i]);
        }
    }
    else
    {
        for (int i = 0; i < np * n - vec_size; i++)
        {
            int val = rand() % 10;
            vec1.push_back(val);
         
            val = rand() % 10;
            vec2.push_back(val);
         
            dot_prod += (vec1[i] * vec2[i]);
        }
    }
 
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
            dot_prod += temp;
            //cout << "pid = " << pid << " dot_prod = " << dot_prod << endl;
        }
        else if (pid % (int)pow(2,j-1) == 0)
        {
            //The sending process sends its data and exits
            MPI_Send (&dot_prod, 1, MPI_INT, pid-j, 0, MPI_COMM_WORLD);
            //cout << pid << " have sent!" << endl;
            break;
        }
        j++;
    }
 
    if (pid == 0)
      cout << "The dot product is " << dot_prod << endl;
 
    MPI_Finalize();

    return 0;
}
