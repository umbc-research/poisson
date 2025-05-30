#pragma once
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
void write_sol_vector(long l_N, long N, double *v, int np, int id, char *filename);
void write_ser_vector(long count, double *v, char *filename);
void writeOut(FILE *outfile, long count, double *v);
