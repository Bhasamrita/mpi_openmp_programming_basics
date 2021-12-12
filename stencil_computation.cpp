#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <float.h>
#include <limits.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>

using namespace std;
using namespace std::chrono;

double difference_of(double new_val, double old_val);

int main(int argc, char** argv)
{
    srand((unsigned) time(0));
 
    MPI_Status status;
 
    //Initialising the MPI environment
    MPI_Init(&argc, &argv);
 
    //Initialing the total number of processes
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
 
    //Assigning a rank to each of the processes
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
 
    int threads, n, ind_n;
    double thres;
    if (pid == 0)
    {
        //Initialising the total number of threads used by each of the processes
        int i = rand() % 6;
        threads = pow(2, i);
        cout << "The number of threads to be used by each process: " << threads << endl;
     
        //Initialising the order of the matrix. The order of the matrix is considered to be varying between 4 and 10000.
        n = rand() % 9997 + 4;
        cout << "The order of the matrix is: " << n << endl;
     
        //Initialising the threshold value or the stopping condition. It ranges from 0 to 1.
        thres = (double) rand() / RAND_MAX;
        cout << "The threshold value is: " << thres << endl;
     
        //Initialising the number of rows computed by individual processes
        ind_n = ceil((double) n / np);
        cout << "The individual rows to be computed by each process (except the last one) is: " << ind_n << endl;
    }

    /*Each process gets the value for the number of threads it needs to launch,
      the order of the matrix (the number of columns each process will hold), the
      threshold value, and the number of rows each process will hold. */
    MPI_Bcast(&threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thres, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ind_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    //For the last process the number of rows it needs to hold may be a little different from the others
    if (pid == np -1)
        ind_n = n - pid * ind_n;
 
    /*The matrix is independently initialised in each of the processes. This part can be visualised as
      each process collecting its part of the matrix that is stored in some storage space. */
    double val;
    vector<double> vec;
    vector<vector<double>> mat;
    
    #pragma parallel for num_threads(threads)
    for (int i = 0; i < ind_n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            val = (double) rand() / RAND_MAX * 100;
            vec.push_back(val);
        }
        mat.push_back(vec);
    }
 
    //Declaring two temporary arrays for each process to hold one row each from the top and bottom.
    //The first and last processes need just one temporary array. Sending arrays are also declared.
    int send_temp1[n], temp_arr1[n], send_temp2[n], temp_arr2[n];
 
    //Except process 0, all processes will be sending their first row to the process just before them.
    //Except process 0, all processes will be receiving the last row of the process just before them.
    if (pid != 0)
    {
        for (int i = 0; i < n; i++)
            send_temp1[i] = mat[0][i];
     
        MPI_Send(&send_temp1, n, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
        MPI_Recv(&temp_arr2, n, MPI_DOUBLE, pid-1, 1, MPI_COMM_WORLD, &status);
    }

    //Except process np-1, all processes will be sending their last row to the process just after them.
    //Except process np-1, all processes will be receiving the first row of the process just after them.
    if (pid != n-1)
    {
        for (int i = 0; i < n; i++)
            send_temp2[i] = mat[ind_n - 1][i];
     
        MPI_Send(&send_temp2, n, MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD);
        MPI_Recv(&temp_arr1, n, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &status);
    }
 
    //Now that all the processes are equipped with all the details that they need, the computation begins.
 
    double max_diff = DBL_MAX, diff = DBL_MAX;
    double left, right, up = 1, down = 1;
    vector<vector<double>> new_val_mat (ind_n, vector<double> (n, -1));
    int iter = 0;
 
    auto start_time = high_resolution_clock::now();
 
    while (max_diff >= thres)
    {
        if (pid == 0) iter++;
     
        #pragma parallel for num_threads(threads)
        for (int i = 0; i < ind_n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == 0 && pid == 0) 
                    up = 0;
                else if (i == 0) 
                    up = temp_arr1[j];
                else 
                    up = temp_arr1[j];
             
                if (i == ind_n-1 && pid == np-1) 
                    down = 0;
                else if (i == ind_n-1) 
                    down = temp_arr2[j];
                else 
                    down = mat[i+1][j];

                if (j == 0) 
                    left = 0;
                else 
                    left = mat[i][j-1];
             
                if (j == n-1) 
                    right = 0;
                else 
                    right = mat[i][j+1];
             
                int new_val = (mat[i][j] + up + down + left + right) / 5;
                diff = difference_of(new_val, mat[i][j]);

                if (max_diff == DBL_MAX)  
                    max_diff = diff;
                if (max_diff < diff)  
                    max_diff = diff;
             
                new_val_mat[i][j] = new_val;
            }
            mat = new_val_mat;
        }
    }
 
    auto ending_time = high_resolution_clock::now();
 
    /*At this point I feel like I could have done better with the above piece of code.
      Something like parallelising it even more, and probably using all the old values at
      the same time and updating the new values at the same time so that the matrix
      new_val_mat isn't required anymore. But I couldn't think of a way to do it till now, so
      I kept this piece of code itself. */
 
    if (pid == 0)
    {
        auto time = duration_cast<microseconds>(ending_time - start_time);
     
        cout << "The time taken for stencil computation is: " << time.count() << endl;
        cout << "The number of iterations taken for stencil computation is: " << iter << endl;
        cout << "The total number of threads launched is: " << threads * np << endl;
    }
 
    //Cleaning up all MPI states before exiting
    MPI_Finalize();
 
    return 0;
}

double difference_of (double new_val, double old_val)
{
    double diff = new_val - old_val;
 
    if (new_val > old_val)
        return diff;
    else
        return -diff;
}
