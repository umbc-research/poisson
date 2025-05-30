#include "main.h"

#define SQ(X) ((X)*(X))

int main (int argc, char *argv[]) {

  int id, np, namelen, idleft, idright;
  int rem, lenm;
  char name[MPI_MAX_PROCESSOR_NAME];
  char filename[64];
  char message[100];
  int iter, maxit, flag;
  long N, l_N, n, l_n, i, j, l_j, l_k, l_m;
  long l_ia, l_ib;
  long l_Ntmp, l_ntmp;
  double *l_u, *l_r, *l_p, *l_q;
  double *x, *y;
  double *gl, *gr;
  double h, tol, relres;
  double enorminf, l_enorminf, err_ij;
  double start, end;
  FILE *idfil;
  MPI_Status status;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Get_processor_name (name, &namelen);

#ifdef BUILD_DATE
  if (id == 0) {
    printf("BUILD DATE: %s\n", BUILD_DATE);
    fflush(stdout);
  }
#endif
#ifdef BUILD_DIR
  if (id == 0) {
    printf("BUILD DIR: %s\n", BUILD_DIR);
    fflush(stdout);
  }
#endif

  if (id > 0)    idleft  = id - 1;
  else           idleft  = MPI_PROC_NULL;
  if (id < np-1) idright = id + 1;
  else           idright = MPI_PROC_NULL;

  /* process command-line inputs: */
  if (argc != 4){
    if (id == 0){
      printf ("Usage: \"poisson N tol maxit\" \n");
      printf (" with integer N, real tol, and integer maxit\n");
      fflush(stdout);
    }
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  N     = (long)(atof(argv[1]));
  tol   =       (atof(argv[2]));
  maxit =  (int)(atof(argv[3]));

  h = 1.0/(N + 1.0);  /*step size*/
  n = N*N; /* = N^2 */

  /* in the following distribution, the first rem=N%np processes get
   * N/np+1 values, and the remaining np-rem ones get N/np,
   * for a total of rem*(N/np+1) + (np-rem)*(N/np)
   * = rem*(N/np) + rem + np*(N/np) - rem*(N/np)
   * = (rem+np-rem)*(N/np) + rem = np*(N/np) + rem = (N-rem) + rem = N */
  rem = N%np; /* remainder when dividing N values into np MPI processes */
  if (id < rem) { /* if id is 0, 1, ..., rem-1 then one more value: */
    l_N = N/np + 1;                  /* local number of values */
  } else {
    l_N = N/np;                      /* local number of values */
  }
  l_n = N * l_N;    /* size of local matrix l_u */
  /*
  if (rem != 0) {
    if (id == 0){
      printf ("Error: np must divide N: N = %5ld, np = %4ld\n", N, np);
      fflush(stdout);
    }
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  */

  /* Notice carefully: start and ending index are coded "C-style",
   * i.e., a for-loop should read: for(l_i=l_ia; l_i<l_ib; l_i++)
   * with a strictly less than comparison. */
  rem = N%np; /* remainder when dividing N values into np MPI processes */
  if (id < rem) { /* if id is 0, 1, ..., rem-1 then one more value: */
    l_ia = id*l_N;                   /* local starting index */
    l_ib = l_ia + l_N;               /* local ending index */
  } else {
    l_ia = rem*(l_N+1)+(id-rem)*l_N; /* local starting index */
    l_ib = l_ia + l_N;               /* local ending index */
  }

  /* test output: */
  if (np <= 8) { /* for more than np=8 this becomes unreadable */
    sprintf(message, "P%03d: l_n=%12ld; l_N=%6ld from l_ia=%6ld to l_ib=%6ld",
                     id, l_n, l_N, l_ia, l_ib);
    if (id == 0) {
      for (i = 0; i < np; i++) {
        if (i > 0)
          MPI_Recv(message,100,MPI_CHAR,i,i,MPI_COMM_WORLD,&status);
        printf("[%3d] %s\n", id, message);
      }
    } else {
      lenm = 1 + strlen(message);
      MPI_Send(message,lenm,MPI_CHAR,0,id,MPI_COMM_WORLD);
    }
  }

  /* allocation and initialization to zero vectors of the 4 large vectors: */
  l_u = allocate_double_vector_openmp(l_n); /* solution vector u */
  l_r = allocate_double_vector_openmp(l_n); /* res. r = right-side b */
  l_p = allocate_double_vector_openmp(l_n); /* search direction p */
  l_q = allocate_double_vector_openmp(l_n); /* auxiliary vector q */

  x = allocate_double_vector_openmp(N);
  y = allocate_double_vector_openmp(N);

  gl = allocate_double_vector_openmp(n/N);
  gr = allocate_double_vector_openmp(n/N);

  for (i = 0; i < N; i++) { x[i] = h * (double)(i+1); }
  for (j = 0; j < N; j++) { y[j] = h * (double)(j+1); }

  for (i = 0; i < n/N; i++) { gl[i] = 0.0; }
  for (i = 0; i < n/N; i++) { gr[i] = 0.0; }

  /* initialization of l_r to values of right-hand side vector b: */
  #ifdef _OPENMP
  #pragma omp parallel for private(i,j,l_k)
  #endif
  for (l_j = 0; l_j < l_N; l_j++){
    j = l_j + l_ia; /* l_j + l_N*id; */
    for (i = 0; i < N; i++){
      l_k = i + N*l_j;
      l_r[l_k] = SQ(h)*(-2*SQ(M_PI))*
                 ( (cos(2*M_PI*x[i]))*SQ(sin(M_PI*y[j]))
                 + SQ(sin(M_PI*x[i]))*(cos(2*M_PI*y[j])) );
    }
  }

  if (id == 0) {
    printf("N = %6ld, tol = %10.1g, maxit = %d\n", N, tol, maxit);
    fflush(stdout);
  }
#ifdef _OPENMP
    #pragma omp parallel
    {
    if (omp_get_thread_num()==0) {
      if (id == 0) printf("Using OpenMP: np = %4d, num_threads = %4d\n",
                          np, omp_get_num_threads());
      }
    }
#else
    if (id == 0) printf("Not using OpenMP: np = %4d\n", np);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  if (id == 0) {
    printf("N = %12ld, np = %4d, l_N = %12ld\n", N, np, l_N);
    printf("n = %12ld, np = %4d, l_n = %12ld\n", n, np, l_n);
#ifdef BLAS
    printf("Using BLAS\n");
#else
    printf("Not using BLAS\n");
#endif
    fflush(stdout);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();  /* start time */

  cg(l_u, &flag, &relres, &iter, l_r, tol, maxit,
     l_p, l_q, l_n, l_N, N, id, idleft, idright, np, MPI_COMM_WORLD, gl, gr);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();  /* end time */

  /*compute L^inf norm of error*/
  l_enorminf = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for private(i,j,l_k)
  #endif
  for (l_j = 0; l_j < l_N; l_j++){
    j = l_j + l_ia; /* l_j + l_N*id; */
    for (i = 0; i < N; i++) {
      l_k = i + N*l_j;
      err_ij = fabs(l_u[l_k] - (SQ(sin(M_PI*x[i]))*SQ(sin(M_PI*y[j]))));
      if (err_ij > l_enorminf)
        l_enorminf = err_ij;
    }
  }
  MPI_Reduce(&l_enorminf,&enorminf,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  if (N <= 5) {
    if (id == 0) {
      for (i = 0; i < np; i++) {
        rem = N%np; /* remainder when dividing N values into np MPI processes */
        if (i < rem) { /* if id is 0, 1, ..., rem-1 then one more value: */
          l_Ntmp = N/np + 1;                  /* local number of values */
        } else {
          l_Ntmp = N/np;                      /* local number of values */
        }
        l_ntmp = N * l_Ntmp;

        if (i == 0) {
          for (l_m=0; l_m<l_ntmp; l_m++)
            l_q[l_m] = l_u[l_m];
        } else {
          MPI_Recv(l_q,l_ntmp,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
        }

        for (l_m=0; l_m<l_ntmp; l_m++)
          printf("%24.16e\n", l_q[l_m]);
      }
    } else { /* id > 0 */
      MPI_Send(l_u,l_n,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }
    fflush(stdout);
  }

  if (id == 0) {
    printf("flag = %1d, iter = %d\n", flag, iter);
    printf("relres             = %24.16e\n", relres);
    printf("h                  = %24.16e\n", h);
    printf("h^2                = %24.16e\n", h*h);
    printf("enorminf           = %24.16e\n", enorminf);
    printf("C = enorminf / h^2 = %24.16e\n", (enorminf/(h*h)));
    printf("wall clock time    = %10.2f seconds\n", end-start);
    diag_time(end-start);
    fflush(stdout);
  }

  /* log nodes used and memory usage */
  nodesused();
  diag_memory();

  if (N < 69) {
    write_sol_vector(l_N, N, l_u, np, id, "solu0000.dat");
    write_ser_vector(N, x, "x.dat");
    write_ser_vector(N, y, "y.dat");
  }

  MPI_Finalize();

  free_vector(l_q);
  free_vector(l_r);
  free_vector(l_p);
  free_vector(l_u);
  free_vector(x);
  free_vector(y);
  free_vector(gl);
  free_vector(gr);

  return 0;
}
