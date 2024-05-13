/*
 * sor.c:
 * Successive over relaxation
 * MPI version implementing a red-black SOR, based on earlier Orca source.
 *
 * Thilo Kielmann, 12/08/1998
 *
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"

int even (int i){
  return ( ( ( i / 2 ) * 2 ) == i ) ? 1 : 0;
}

double abs_d (double d){
  return (d < 0.0) ? -d : d;
}

double stencil (double**G, int row, int col){
  return ( G[row-1][col] + G[row+1][col] + G[row][col-1] + G[row][col+1] )
         / 4.0;
}

/* compute lower and upper bound of the grid stripe */
static void
get_bounds(int *lb, int *ub, int n, int size, int rank)
{
    int nlarge, nsmall;
    int size_large, size_small;

    nlarge = n % size;
    nsmall = size - nlarge;

    size_small = n / size;
    size_large = size_small + 1;

    if (rank < nlarge) {            /* I'll  have a large chunk */
	*lb = rank * size_large;
	*ub = *lb + size_large;
    } else {
	*lb = nlarge * size_large + (rank - nlarge) * size_small;
	*ub = *lb + size_small;
    }
}

int main (int argc, char *argv[]){

  int N;                            /* problem size */
  int ncol,nrow;                    /* number of rows and columns */
  double TOLERANCE = 0.0002;        /* termination criterion */
  double Gnew,r,omega,              /* some float values */
         stopdiff,maxdiff,diff;     /* differences btw grid points in iters */
  double **G;                       /* the grid */
  int i,j,phase,iteration;          /* counters */
  double t_start,t_end;             /* time values */
  int size,rank,pred,succ;          /* process ranks */
  int lb,ub;                        /* lower and upper bound of grid stripe */
  MPI_Status status;                /* dummy for MPI_Recv */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /* ranks of predecessor and successor for row exchanges */
  pred = (rank == 0) ? MPI_PROC_NULL : rank - 1;
  succ = (rank == size-1) ? MPI_PROC_NULL : rank + 1;

  /* set up problem size */
  N = 0;
  if ( argc > 1 )
    N = atoi(argv[1]);
  if (N == 0)
    N = 1000;

  if ( N < size ) N = size; /* give each process at least one row */

  if (rank == 0 ) {
      printf("Running SOR on %d nodes with %d rows\n", size, N);
  }

  N += 2; /* add the two border lines */
          /* finally N*N (from argv) array points will be computed */

  /* set up a quadratic grid */
  ncol = nrow = N;
  r        = 0.5 * ( cos( M_PI / ncol ) + cos( M_PI / nrow ) );
  omega    = 2.0 / ( 1 + sqrt( 1 - r * r ) );
  stopdiff = TOLERANCE / ( 2.0 - omega );
  omega   *= 0.8;                   /* magic factor */

  /* get my stripe bounds and malloc the grid accordingly */
  get_bounds(&lb, &ub, N-1, size, rank);
  if (lb == 0) lb = 1; /* row 0 is static */

  G = (double**)malloc(nrow*sizeof(double*));
  if ( G == 0 ){ printf("malloc failed\n"); exit(42); }
  for (i = lb-1; i<=ub; i++){ /* malloc the own range plus one more line */
                              /* of overlap on each border */
    G[i] = (double*)malloc(ncol*sizeof(double));
    if ( G[i] == 0 ){ printf("malloc failed\n"); exit(42); }
  }

  /* initialize the grid */
  for (i = lb-1; i <= ub; i++){
    for (j = 0; j < ncol; j++){
      if (i == 0) G[i][j] = 4.56;
      else if (i == nrow-1) G[i][j] = 9.85;
      else if (j == 0) G[i][j] = 7.32;
      else if (j == ncol-1) G[i][j] = 6.88;
      else G[i][j] = 0.0;
    }
  }

  /* now do the "real" computation */
  t_start = MPI_Wtime();
  iteration = 0;
  do {
    MPI_Send(G[lb]  ,ncol,MPI_DOUBLE,pred,42,MPI_COMM_WORLD);
    MPI_Send(G[ub-1],ncol,MPI_DOUBLE,succ,42,MPI_COMM_WORLD);
    MPI_Recv(G[lb-1],ncol,MPI_DOUBLE,pred,42,MPI_COMM_WORLD,&status);
    MPI_Recv(G[ub]  ,ncol,MPI_DOUBLE,succ,42,MPI_COMM_WORLD,&status);
    maxdiff = 0.0;
    for ( phase = 0; phase < 2 ; phase++){
      for ( i = lb ; i < ub ; i++ ){
        for ( j = 1 + (even(i) ^ phase); j < ncol-1 ; j += 2 ){
          Gnew = stencil(G,i,j);
          diff = abs_d(Gnew - G[i][j]);
          if ( diff > maxdiff )
            maxdiff = diff;
          G[i][j] = G[i][j] + omega * (Gnew-G[i][j]);
        }
      }
    }
    diff = maxdiff;
    MPI_Allreduce(&diff,&maxdiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    iteration++;
  } while (maxdiff > stopdiff);
  t_end = MPI_Wtime();

  if ( rank == 0 ){
    printf("SOR %d x %d took %f s\n",nrow-2,ncol-2,t_end-t_start);
    printf("using %d iterations, diff is %f (allowed diff %f)\n",
           iteration,maxdiff,stopdiff);
  }

  MPI_Finalize();
  return 0;
}
