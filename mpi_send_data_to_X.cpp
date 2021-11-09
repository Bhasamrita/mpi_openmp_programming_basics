/*The program attempts to exchange array elements.
  First, process 0 becomes the sending process and generates random number / np
  to find the receiving process. Then it sends its array elements to the receiving
  process. The receiving process adds these incoming elements into its array.
  For the next iteration the receiving process becomes the sender. This is 
  continued for a predetermined number of iterations. The array elements contained
  in all the processes are displayed in the end. */

#include <iostream>
#include <ctime>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    srand((unsigned) time(0));
 
    MPI_Status status;
 
    MPI_Init(&argc, &argv);
 
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
 
    int iterations, arr_size;
    if (pid == 0)
    {
        //Initialising the array size for each process
        arr_size = rand() % 6 + 5;
        cout << "The array size for each process is " << arr_size << endl;
        //Initialising the number of iterations
        iterations = rand() % 6 + 5;
        cout << "The number of iterations is " << iterations << endl << endl;
    }

    //Broadcasting the array size to all processes from process 0
    MPI_Bcast(&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //Broadcasting the number of iterations to all processes from process 0
    MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    //Every process then initialises a vector of size arr_size
    vector<int> vec;
    int val;
    for (int i = 0; i < arr_size; i++)
    {
        val = rand() % 10;
        vec.push_back(val);
        //cout << vec[i] << " ";
    }
    //cout << endl;
 
    MPI_Barrier(MPI_COMM_WORLD);
 
    //The exchange of array elements starts here
    int send_pid = 0, recv_pid, send_size, recv_size;
    while(iterations--)
    {
        if (pid == send_pid)
            recv_pid = rand() % np;
     
        MPI_Bcast(&recv_pid, 1, MPI_INT, send_pid, MPI_COMM_WORLD);
     
        if (pid == send_pid)
        {
            send_size = (int)vec.size();
            int temp[send_size];
            for (int i = 0; i < vec.size(); i++)
              temp[i] = vec[i];
         
            MPI_Send(&send_size, 1, MPI_INT, recv_pid, 0, MPI_COMM_WORLD);
            MPI_Send(&temp, send_size, MPI_INT, recv_pid, 1, MPI_COMM_WORLD);
         
            //cout << "send_pid = " << send_pid << " and send_size = " << send_size << endl << endl;
        }
        
        if (pid == recv_pid)
        {
            MPI_Recv(&recv_size, 1, MPI_INT, send_pid, 0, MPI_COMM_WORLD, &status);
         
            int temp[recv_size] = {0};
            MPI_Recv(&temp, recv_size, MPI_INT, send_pid, 1, MPI_COMM_WORLD, &status);
         
            //cout << "recv_pid = " << recv_pid << " and recv_size = " << recv_size << endl << endl;
         
            for (int i = 0; i < recv_size; i++)
            {
                vec.push_back(temp[i]);
                //cout << temp[i] << " ";
            }
            //cout << endl << endl;
            send_pid = recv_pid;
        }
     
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&send_pid, 1, MPI_INT, send_pid, MPI_COMM_WORLD);
    }
 
    if (pid == 0)
      cout << "The elements in all the processes are - " << endl;
    MPI_Barrier(MPI_COMM_WORLD);
 
    for (auto x : vec)
      cout << x << " ";
    cout << endl;
 
    MPI_Finalize();
 
    return 0;
}
