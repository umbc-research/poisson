#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "utilities.h"
#include "cg.h"
#include "check_memory.h"
#include "nodesused.h"
#include "diag_time.h"
#include "outputs.h"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BLAS
#include <mkl.h>
#endif

#ifndef M_PI
#define M_PI          (3.1415926535897932384626433832795029)
#endif

#endif

