#ifndef MEMORY_H
#define MEMORY_H

#include <stdio.h>
#include <stdlib.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

int *allocate_int_vector (int n);
long int *allocate_long_int_vector (long int n);
double *allocate_double_vector (long n);
double *allocate_double_vector_openmp (long n);
double *allocate_double_vector_openmp_nz (long n, double val);
long *allocate_long_int_vector_openmp (long n);
void free_vector (void *x);


#endif

