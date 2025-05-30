/* This file prints the output for the poisson code is the mesh 
 * is sufficiently small.
 *
 */
#include "outputs.h"
void write_sol_vector(long l_N, long N, double *v, int np, int id, char *filename) {
    /*char *tStep = "0000";*/
    /*char *filename = sprintf("solu%s.dat", tStep);*/
    long n = l_N * N;
    // Allocate the RECV array
    FILE *outfile = fopen(filename, "w");
#ifdef PARALLEL
    MPI_Status status;
    // Process 0 always has  the largest buffer
    double *rBuffer = allocate_double_vector_openmp(n);
    if (id == 0) {
        // Print out process 0's data
        writeOut(outfile, n, v);
        for (int i = 1; i < np; i++) {
            // Buffer size isn't always the same
            MPI_Recv(&n, 1, MPI_LONG, i, 0, MPI_COMM_WORLD, &status);
            // Fill buffer
            MPI_Recv(rBuffer, n, MPI_DOUBLE, i, 0,  MPI_COMM_WORLD, &status);
            writeOut(outfile, n, rBuffer);
        }
    }
    else {
        // Send data to process 0
        MPI_Send(&n, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(rBuffer, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
#else
    writeOut(outfile, n, v);
#endif
    fclose(outfile);
}

void write_ser_vector(long count, double *v, char *filename) {
    FILE *outfile  = fopen(filename, "w");
    writeOut(outfile, count, v);
    fclose(outfile);
}

void writeOut(FILE *outfile, long count, double *v) {
    long i;
    for (i = 0;  i < count; i++) {
        fprintf(outfile, "%24.16e\n", v[i]);
    }
}


