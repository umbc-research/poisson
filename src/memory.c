#include "memory.h"

/* 09/12/02-10/10/02, updated 02/07/08, 02/23/16 by Matthias K. Gobbert */

int *allocate_int_vector (int n)
{
  int *x;

  x = (int*) malloc (n* sizeof(int));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}


long int *allocate_long_int_vector (long int n)
{
  long int *x;

  x = (long int*) malloc (n* sizeof(long int));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

double *allocate_double_vector (long n)
{
  double *x;

  x = (double*) malloc (n* sizeof(double));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

void free_vector (void *x)
{
  free (x);
}

// Additional allocations to use OpenMP
/*#ifdef _OPENMP*/
double *allocate_double_vector_openmp(long n) {
  double *x = (double *) malloc(n * sizeof(double));
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (long l_k = 0; l_k < n; l_k++) {
    x[l_k] = 0.0;
  }
  return x;
}

// nz = non-zero
double *allocate_double_vector_openmp_nz(long n, double val) {
  double *x = (double *) malloc(n * sizeof(double));
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (long l_k = 0; l_k < n; l_k++) {
    x[l_k] = val;
  }
  return x;
}


long *allocate_long_int_vector_openmp(long n) {
  long *x = (long *) malloc(n * sizeof(long));
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (long l_k = 0; l_k < n; l_k++) {
    x[l_k] = 0;
  }
  return x;
}
/*#endif*/
