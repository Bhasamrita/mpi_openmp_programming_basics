/*Is there life anywhere else in the universe? Are there intelligent civilisations
  anywhere else in the Milky Way Galaxy? These are hyped questions that get asked
  whenever the discussion moves to astronomy and the universe. But the general
  discussions (almost always) lead to nowhere.

  However, in astrophysics, there is a probabilistic equation called Drake's
  equation which gives us the probability that intelligent life exist elsewhere
  in the Milky Way Galaxy -
    N = R_star * fp * ne * fl * fi * fc * L
  where N is the number of civilisations with which humans could communicate
        R_star is the mean rate of star formation
        fp is the fraction of stars that have planets
        ne is the mean number of planets that could support life per star with
           planets
        fl is the fraction of life supporting planets that develop life
        fi is the fraction of planets with life where life develops intelligence
        fc is the fraction of intelligent civilisations that develops communication
        L is the mean length of time that civilisation can communicate

  The last three factors - fi, fc and L are very hard to predict. This equation,
  although not used to precisely calculate an answer, is heavily used in
  simulation, and probabilistic calculations in astronomy and astrophysics.

  This program is an attempt to parallelize the core part of the problem for faster
  and more effective simulation.
  Formation of the problem (perhaps in an overly simplified manner)-
  1. Let R be the vector which contains the estimated rate of formation of stars
     in each of the nebulae present in the Milky Way Galaxy. It will be a
     20000X1 vector, corresponding to the ~20000 nebulae found in our galaxy.
     This vector contains decimal number between 1 and 2, meaning the mean rate
     of formation of stars in each of the nebulae is a little over 1 star/year.
     This vector is considered to be already given to us.
  2. Let Sp be the binary vector which contains 1 corresponding to stars with
     planets, and 0 corresponding to stars without planets. It will be a
     100000000000X1 (100 billion dimensional) vector, corresponding to the
     100 billion stars in the Milky Way Galaxy. This vector is also considered to
     be given to us.
  3. Let Pl be the vector which contains the number of planets corresponding to
     each star with planets. It is a 40000000000X1 (40 billion dimensional) vector.
     This vector contains integers between 1 and 10, denoting number of planets
     for each star with planets.
  4. Let Lf be the vector that contains the number of life supporting planets per
     star with planets. It is also a 40000000000X1 vector. It contains integers
     between 0 and 5, meaning that there can be a maximum of 5 habitable planets
     for every star with planets.
  5. The vector Dl contains the number of planets that are habitable. It is a
      bool vector that has either 1 or 0 corresponding to every planet having
      life or not. It is a 300000000X1 (300 million dimensional) vector.
  6. Let Li be the vector that contains the number of planets with potential life
     in it. It is a bool vector containing 1 or 0 corresponding to whether a
     particular life supporting planet can host intelligent life or not. The
     length of this vector is not fixed and is estimated.
  7. Let Lc be the vector that contains the number of planets with intelligent
     life. It is bool vector containing 1 or 0 corresponding to whether a
     particular intelligent species can communicate with the universe using radio
     technology within their span of existence or not. The length of this vector
     is estimated.
  8. Let Tc be the vector that contains the estimated length of time that each
     civilisation can communicate. (According to Carl Sagan, the famous cosmologist,
     an intelligent civilisation at one point of time will develop enough
     intelligence to finally bring their own end.) The length is estimated.

  Once the vectors are formed and calculated, the individual terms of this
  equation can be calculated.
  1. R_star will be the average of the values contained in vector R.
  2. fp will be length of Pl vector by length of Sp vector.
  3. ne will be the average of the values contained in the vector Lf.
  4. fl will be length of Li vector by length of Dl vector.
  5. fi will be length of Lc vector by lenght of Li vector.
  6. fc will be length of Tc vector by length of Lc vector.
  7. L will be the average of all the values contained in vector Tc.   */

#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <ctime>
#include <cmath>
#include <vector>

#define R_len 20000
#define Sp_len 100000000000
#define Pl_len 40000000000
#define Lf_len 40000000000
#define Dl_len 3000000000
#define nt 8

using namespace std;

void initialise_bool_vector(vector<bool>& vec, long long len, int np, int pid);
void initialise_int_vector(vector<int>& vec, long long len, int d, int c, int np, int pid);
void initialise_double_vector(vector<double>& vec, long long len, int np, int pid);
void sum_bool_individual (vector<bool>& vec, long long &sum, int np, int pid);
void sum_int_individual (vector<int>& vec, long long &sum, int np, int pid);
void sum_double_individual (vector<double>& vec, double &sum, int np, int pid);
void sum_from_all_processes (long long &sum, int np, int pid);
void double_sum_from_all_processes (double &sum, int np, int pid);
int power(int a, int b);

int main(int argc, char** argv)
{
    //Initialising the rand() function to system time
    srand ((unsigned) time(0));

    //Initialising the MPI environment
    MPI_Init(&argc, &argv);

    //Rank or process id
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Size or number of processes
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //Thread id
    int tid;

    /*Phase 1 :
    Vectors R and Sp are considered to be already initialised and then
    given to us. */
    vector<double> R;
    vector<bool> Sp;
    initialise_double_vector(R, R_len, np, pid);
    initialise_bool_vector(Sp, Sp_len, np, pid);

    /*The vector Sp is computed to find out the number of stars that contain
    planets. That number will be used in phase 2.*/
    long long planets = 0;
    sum_bool_individual(Sp, planets, np, pid);
    sum_from_all_processes (planets, np, pid);
    if (pid == 0)
      cout << "There are " << planets << " stars that contain planets." << endl;

    /*Phase 2:
    After getting the number of stars that contain planets, another vector Pl is
    created. It contains the number of planets corresponding to each star with
    planets.This vector is also initialised before it is given to the program. */

    vector<int> Pl;
    //long long Pl_len = planets;
    initialise_int_vector(Pl, Pl_len, 10, 1, np, pid);

    /*Phase 3:
    Each planet corresponding to each star is experimented to find out if it is
    in the habitable zone. Another vector Lf is created which contains the number
    of planets corresponding to each star with planets, which are in the
    habitable zone. This vector is also initialised before it is given to the
    program. */

    vector<int> Lf;
    //long long Lf_len = planets;
    initialise_int_vector(Lf, Lf_len, 6, 0, np, pid);

    /*Given the vector Lf, the program needs to calculate the total number of
    planets that are in the habitable zone. This number will be used in phase 4. */
    long long habitable = 0;
    sum_int_individual (Lf, habitable, np, pid);
    sum_from_all_processes (habitable, np, pid);
    if (pid == 0)
      cout << "The are a total of " << habitable << " habitable planets." << endl;

    /*Phase 4:
    For every habitable planet, experiments are done to find out if the planet
    can develop life (the results here will not be very accurate). The vector Dl
    contains either 1 or 0 corresponding to whether a planet can develop life
    or not. This vector is initialised before it is given to us. */

    vector<bool> Dl;
    //long long Dl_len = habitable;
    initialise_bool_vector(Dl, Dl_len, np, pid);

    /*Given the vector Dl, the program needs to calculate the number of planets
    that can develop life. This number will be used in phase 5. */
    long long life = 0;
    sum_bool_individual (Dl, life, np, pid);
    sum_from_all_processes (life, np, pid);
    if (pid == 0)
      cout << "A total of " << life << " habitable planets can develop life." << endl;

    /*Phase 5:
    For every planet that can develop life, experiments are done to find out if
    it can develop intelligent life. The vector Li contains this information.
    (The experiments done will involve probabilities and will not be accurate.)
    It is initialised before it is given to the program. */

    vector<bool> Li;
    long long Li_len;
    if (pid == 0)
    {
        Li_len = life;
        MPI_Bcast(&Li_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    }
    else
      MPI_Bcast(&Li_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    //cout << "Li_len = " << Li_len << endl;
    initialise_bool_vector(Li, Li_len, np, pid);

    /*Given the vector Li, the program needs to find out the total number of
    planets that can develop intelligent life. This number will be useful in
    phase 6. */
    long long intelligent = 0;
    sum_bool_individual (Li, intelligent, np, pid);
    sum_from_all_processes(intelligent, np, pid);
    if (pid == 0)
      cout << "There are " << intelligent << " planets that can develop " <<
              "intelligent life." << endl;

    /*Phase 6:
    Given the total number of planets that can develop intelligent life,
    experiments are done to find out which of these planets can develop life
    that is able to communicate to the universe using radio waves. The vector
    Lc contains this information. It is initialised before it is given to the
    program. */

    vector<bool> Lc;
    long long Lc_len;
    if (pid == 0)
    {
        Li_len = intelligent;
        MPI_Bcast(&Lc_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    }
    else
      MPI_Bcast(&Lc_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    //cout << "Lc_len = " << Lc_len << endl;
    initialise_bool_vector(Lc, Lc_len, np, pid);

    /*Given the vector Lc, we need to find out the total number of planets that
    can develop intelligent life that can communicate to the universe using
    radio waves. This number will be useful for phase 7. */
    long long communicate = 0;
    sum_bool_individual (Lc, communicate, np, pid);
    sum_from_all_processes(communicate, np, pid);
    if (pid == 0)
      cout << "There are " << communicate << " planets that can develop " <<
              "intelligent life that can communicate to the universe " <<
              "using radio waves." << endl;

    /*Phase 7:
    Given the number of planets that can develop intelligent life that can
    communicate, experiments are done to find out how many years can they
    communicate before they bring their own doom. (The experiments done contain
    probabilistic values and are not accurate.) The vector Tc is initiased with
    this values before it is given to the program. */

    vector<int> Tc;
    long long Tc_len;
    if (pid == 0)
    {
        Tc_len = communicate;
        MPI_Bcast(&Tc_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    }
    else
      MPI_Bcast(&Tc_len, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    //cout << "Tc_len = " << communicate << endl;
    initialise_int_vector(Tc, Tc_len, 400, 101, np, pid);

    /*Equipped with every necessay information, the program can now calculate
    the individual terms present in Drake's equation. */

    /*Calculating the term R_star:
    The vector R contains the mean values of star formation for each nebulae in
    the Milky Way. So, R_star will be the average of these values. */

    double R_star;
    sum_double_individual(R, R_star, np, pid);
    double_sum_from_all_processes(R_star, np, pid);
    if (pid == 0)
    {
        R_star /= R_len;
        cout << "R_star = " << R_star << endl;
    }

    /*Calculating the term fp:
    The vector Sp contains information on every star in the Milky Way. Whether
    each star contains planet or not. So, fp will be summation of values in
    Sp divided by length of Sp. */

    double fp;
    if (pid == 0)
    {
        fp = (double) planets / Sp_len;
        cout << "fp = " << fp << endl;
    }

    /*Calculating the term ne:
    The vector Lf contains the number of habitable planets corresponding to each
    star with planets. So, ne will be the average of all the values contained in
    the vector Lf. */

    double ne;
    if (pid == 0)
    {
        ne = (double) habitable / Lf_len;
        cout << "ne = " << ne << endl;
    }

    /*Calculating the term fl:
    Vector Dl contains the number of habitable planets that can develop life. So,
    fl will be total number of 1s in Dl divided by the length of Dl. */

    double fl;
    if (pid == 0)
    {
        fl = (double) life / Dl_len;
        cout << "fl = " << fl << endl;
    }

    /*Calculating the term fi:
    The vector Li contains the number of life-supporting planets that can develop
    intelligent life. So, fi will be the number of 1s in Li divided by the length
    of Li. */

    double fi;
    if (pid == 0)
    {
        fi = (double) intelligent / Li_len;
        cout << "fi = " << fi << endl;
    }

    /*Calculating the term fc:
    The vector Lc contains the number of intelligent life-supporting planets that
    can develop communication. So, fc will be the number of 1s in Lc divided by
    the length of Lc. */

    double fc;
    if (pid == 0)
    {
        fc = (double) communicate / Lc_len;
        cout << "fc = " << fc << endl;
    }

    /*Calculating the term L:
    The vector Tc contains the mean number of years that each planet with
    communication abled intelligent life can communicate for. So, L will be the
    average of all these values. */

    double L;
    long long sum = 0;
    sum_int_individual (Tc, sum, np, pid);
    sum_from_all_processes(sum, np, pid);
    if (pid == 0)
    {
        L = (double) sum / Tc_len;
        cout << "time = " << L << endl;
    }

    /*Finally, calculating Drake's equation:
    Equipped with all the values for the individual terms involved in this
    equation, the program can now calculate for N. */

    if (pid == 0)
    {
        double N;
        N = (double) R_star * fp * ne * fl * fi * fc * L;
        cout << "N = " << N << endl;
    }

    //Cleaning up all MPI states before exiting
    MPI_Finalize();

    return 0;
}

/*All the function definitions start from here. */

void initialise_bool_vector (vector<bool>& vec, long long len, int np, int pid)
{
    long long n = ceil(len / np);
    int val;
    if (pid != np - 1)
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < n; i++)
        {
            val = rand() % 2;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
    else
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < len - pid*n; i++)
        {
            val = rand() % 2;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
}

void initialise_int_vector (vector<int>& vec, long long len, int d, int c, int np, int pid)
{
    long long n = ceil(len / np);
    int val;
    if (pid != np - 1)
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < n; i++)
        {
            val = rand() % d + c;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
    else
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < len - pid*n; i++)
        {
            val = rand() % d + c;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
}

void initialise_double_vector (vector<double>& vec, long long len, int np, int pid)
{
    long long n = ceil(len / np);
    double val;
    if (pid != np - 1)
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < n; i++)
        {
            val = (double) rand() / RAND_MAX + 1;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
    else
    {
        #pragma omp parallel for num_threads(nt)
        for (long long i = 0; i < len - pid*n; i++)
        {
            val = (double) rand() / RAND_MAX + 1;
            vec.push_back(val);
            //cout << val << " ";
        }
        //cout << endl;
    }
}

void sum_bool_individual (vector<bool>& vec, long long &sum, int np, int pid)
{
    #pragma omp parallel for num_threads(nt)
    for (auto i : vec)
      sum += i;
}

void sum_int_individual (vector<int>& vec, long long &sum, int np, int pid)
{
    #pragma omp parallel for num_threads(nt)
    for (auto i : vec)
      sum += i;
}

void sum_double_individual (vector<double>& vec, double &sum, int np, int pid)
{
    #pragma omp parallel for num_threads(nt)
    for (auto i : vec)
      sum += i;
}

void sum_from_all_processes (long long &sum, int np, int pid)
{
    MPI_Status status;

    /*If we have 8 processes, then the below code ensures the following -
      - In the first iteration, all even processes will receive,
        and odd processes will send and exit the loop.
      - In the second iteration, all 4modulo processes will receive,
        and the remaining will send and exit.
      - In the third iteration, all 8modulo processes will receive,
        and the remaining will send and exit.
      This loop will run for log_base_2(np) times.
      The objective here is to add up the numbers obtained by all the processes
      so that we get the number of stars that have planets. */

    int j = 1;
    while (power(2, j) <= np)
    {
        //cout << "pid = " << pid << " j = " << j << endl;
        if (pid % power(2,j) == 0)
        {
            /*The receiving process initialises a buffer to receive the data
              and then adds it to its local sum variable. */
            long long temp;
            MPI_Recv (&temp, 1, MPI_LONG_LONG, pid+j, 0, MPI_COMM_WORLD, &status);
            sum += temp;
            //cout << "pid = " << pid << " local_sum = " << local_sum << endl;
        }
        else if (pid % power(2,j-1) == 0)
        {
            MPI_Send (&sum, 1, MPI_LONG_LONG, pid-j, 0, MPI_COMM_WORLD);
            //cout << pid << " have sent!" << endl;
            break;
        }
        j++;
    }
}

void double_sum_from_all_processes (double &sum, int np, int pid)
{
    MPI_Status status;
    int j = 1;
    while (power(2, j) <= np)
    {
        //cout << "pid = " << pid << " j = " << j << endl;
        if (pid % power(2,j) == 0)
        {
            /*The receiving process initialises a buffer to receive the data
              and then adds it to its local sum variable. */
            double temp;
            MPI_Recv (&temp, 1, MPI_DOUBLE, pid+j, 0, MPI_COMM_WORLD, &status);
            sum += temp;
            //cout << "pid = " << pid << " local_sum = " << local_sum << endl;
        }
        else if (pid % power(2,j-1) == 0)
        {
            MPI_Send (&sum, 1, MPI_DOUBLE, pid-j, 0, MPI_COMM_WORLD);
            //cout << pid << " have sent!" << endl;
            break;
        }
        j++;
    }
}

int power(int a, int b)   //Defining the power() function that returns int
{
    int p = 1;
    for (int i = 0; i < b; i++)
      p *= a;
    return p;
}

//P.S. The program was tested for a maximum vector length of 1000000 (1 million)

