/*
 *
 *  Optimization of the Sieve of Eratosthenes:
 *  a) not consider even numbers (as they are all multiply of 2 and this will save 50% memory usage)
 *  b) all the process compute the sieve (this removes a broadcast communication)
 *
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[]) {
    long count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    long first;        /* Index of first multiple */
    long global_count; /* Global prime count */
    long high_value;   /* Highest value on this proc */
    int i;
    int rank;           /* Process ID number */
    long index;        /* Index of current prime */
    long low_value;    /* Lowest value on this proc */
    char *marked;       /* Portion of 3,...,'n' */
    long n;            /* Sieving from 3, ..., 'n' */
    int processes;            /* Number of processes */
    long proc0_size;   /* Size of proc 0's subarray */
    long prime;        /* Current prime */
    long size;         /* Elements in 'marked' */

    MPI_Init(&argc, &argv);


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);


    /* Start the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2) {
        if (!rank){
            fprintf (stdout,"./%s <m>\n", argv[0]);
            fflush(stdout);
        }
        MPI_Finalize();
        exit (1);
    }

    n = strtol(argv[1], NULL, 10);

    low_value = (long) 2 + rank*(n-1)/processes;
    high_value = (long) 1 + (rank+1)*(n-1)/processes;
    size = high_value - low_value + 1;

    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = (long) (n-1)/processes;

    if ((2 + proc0_size) < (long) sqrt((double) n)) {
        if (!rank){
            fprintf (stdout,"Too many processes\n");
        }
        MPI_Finalize();
        exit (1);
    }

    marked = (char *) malloc(size);
    if (!marked) {
        fprintf (stdout,"Cannot allocate enough memory\n");
        MPI_Finalize();
        exit (1);
    }

    for (i = 0; i < size; i++)
        marked[i] = 0;

    if (!rank)
        index = 0;

    prime = 3;

    do {
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else {
            if (!(low_value % prime))
                first = 0;
            else
                first = prime - (low_value % prime);
        }

        for (i = first; i < size; i += prime)
            marked[i] = 1;

            while (marked[++index])
                ;
            prime = index + 2;
        //fprintf(stdout,"p %d prime: %d\n",rank,prime);

    } while (prime * prime <= n);

    count = 0;

    for (i = 0; i < size; i++)
        if (!marked[i])
            count++;

    if (processes > 1) MPI_Reduce (&count, &global_count, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


    /* Stop the timer */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();


    /* Print the results */

    if (!rank) {
        fprintf (stdout, "There are %ld primes less than or equal to %ld\n", global_count, n);
        fprintf (stdout, "SIEVE (%d) %10.6f\n", processes, elapsed_time);
    }

    MPI_Finalize ();
    return 0;
}

