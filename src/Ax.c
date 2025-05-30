#include "Ax.h"

void Ax(double *l_v, double *l_u, long l_n, long l_N, long N,
        int id, int idleft, int idright, int np, MPI_Comm comm,
        double *gl, double *gr) {

  long i, l_j;
  double tmp;
  MPI_Status statuses[4];
  MPI_Request req[4];

  if (np > 1) {
    MPI_Irecv(gl             ,N,MPI_DOUBLE,idleft ,1,MPI_COMM_WORLD,&(req[0]));
    MPI_Irecv(gr             ,N,MPI_DOUBLE,idright,2,MPI_COMM_WORLD,&(req[1]));
    MPI_Isend(&l_u[N*(l_N-1)],N,MPI_DOUBLE,idright,1,MPI_COMM_WORLD,&(req[2]));
    MPI_Isend(&l_u[0]        ,N,MPI_DOUBLE,idleft ,2,MPI_COMM_WORLD,&(req[3]));
  }

  #ifdef _OPENMP
  #pragma omp parallel for private(i,tmp)
  #endif
  for (l_j = 1; l_j < l_N-1; l_j++) {
    for (i = 0; i < N; i++){
                tmp = 4.0*l_u[i   + N* l_j   ];
                   tmp -= l_u[i   + N*(l_j-1)];
      if (i >  0 ) tmp -= l_u[i-1 + N* l_j   ];
      if (i < N-1) tmp -= l_u[i+1 + N* l_j   ];
                   tmp -= l_u[i   + N*(l_j+1)];
      l_v[i+N*l_j] = tmp;
    }
  }

  if (np > 1) {
    MPI_Waitall(4,req,statuses);
  }

  if(l_N == 1){
    l_j = 0;
    #ifdef _OPENMP
    #pragma omp parallel for private(tmp)
    #endif
    for (i = 0; i < N; i++) {
                tmp = 4.0*l_u[i   + N* l_j   ];
      if (id>  0 ) tmp -=  gl[i];
      if (i >  0 ) tmp -= l_u[i-1 + N* l_j   ];
      if (i < N-1) tmp -= l_u[i+1 + N* l_j   ];
      if (id<np-1) tmp -=  gr[i];
      l_v[i+N*l_j] = tmp;
    }
  }
  else{
    l_j = 0;
    #ifdef _OPENMP
    #pragma omp parallel for private(tmp)
    #endif
    for (i = 0; i < N; i++) {
                tmp = 4.0*l_u[i   + N* l_j   ];
      if (id>  0 ) tmp -=  gl[i];
      if (i >  0 ) tmp -= l_u[i-1 + N* l_j   ];
      if (i < N-1) tmp -= l_u[i+1 + N* l_j   ];
                   tmp -= l_u[i   + N*(l_j+1)];
      l_v[i+N*l_j] = tmp;
    }
    l_j = l_N - 1;
    #ifdef _OPENMP
    #pragma omp parallel for private(tmp)
    #endif
    for (i = 0; i < N; i++) {
                tmp = 4.0*l_u[i   + N* l_j   ];
                   tmp -= l_u[i   + N*(l_j-1)];
      if (i >  0 ) tmp -= l_u[i-1 + N* l_j   ];
      if (i < N-1) tmp -= l_u[i+1 + N* l_j   ];
      if (id<np-1) tmp -=  gr[i];
      l_v[i+N*l_j] = tmp;
    }
  }
}
