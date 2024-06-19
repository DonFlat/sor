/*
 * Successive over relaxation
 * (red-black SOR)
 *
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <sys/time.h>

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


int main (int argc, char *argv[]){

  int N;                            /* problem size */
  int ncol,nrow;                    /* number of rows and columns */
  double TOLERANCE = 0.0002;        /* termination criterion */
  double Gnew,r,omega,              /* some float values */
         stopdiff,maxdiff,diff;     /* differences btw grid points in iters */
  double **G;                       /* the grid */
  int i,j,phase,iteration;          /* counters */

  /* set up problem size */
  N = 0;
  if ( argc > 1 )
    N = atoi(argv[1]);
  if (N == 0)
    N = 1000;


  printf("Running SOR with %d rows\n", N);

  N += 2; /* add the two border lines */
          /* finally N*N (from argv) array points will be computed */

  /* set up a quadratic grid */
  ncol = nrow = N;
  r        = 0.5 * ( cos( M_PI / ncol ) + cos( M_PI / nrow ) );
  omega    = 2.0 / ( 1 + sqrt( 1 - r * r ) );
  stopdiff = TOLERANCE / ( 2.0 - omega );
  omega   *= 0.8;                   /* magic factor */

  G = (double**)malloc(N*sizeof(double*));
  if ( G == 0 ){ printf("malloc failed\n"); exit(42); }
  for (i = 0; i<N; i++){ /* malloc the own range plus one more line */
                         /* of overlap on each side */
    G[i] = (double*)malloc(N*sizeof(double));
    if ( G[i] == 0 ){ printf("malloc failed\n"); exit(42); }
  }

  /* initialize the grid */
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      if (i == 0) G[i][j] = 4.56;
      else if (i == N-1) G[i][j] = 9.85;
      else if (j == 0) G[i][j] = 7.32;
      else if (j == N-1) G[i][j] = 6.88;
      else G[i][j] = 0.0;
    }
  }

  struct timeval start, end;
  long mtime, seconds, useconds;
  /* now do the "real" computation */
  iteration = 0;
  gettimeofday(&start, NULL);
  do {
    maxdiff = 0.0;
    // [0, 1]
    for ( phase = 0; phase < 2 ; phase++){
      // [1, second last element]
      for ( i = 1 ; i < N-1 ; i++ ){
        // 1 + (even(1) ^ 1) => 1 + 1 = 2; j => 4, 6, 
        // 1 + (even(1) ^ 0) => 1 + 0 = 1; j => 3
        // Phase 0: j is odd, Phase 1: j is even
        for ( j = 1 + (even(i) ^ phase); j < N-1 ; j += 2 ){
          Gnew = stencil(G,i,j); // average of four direction
          diff = abs_d(Gnew - G[i][j]);
          if ( diff > maxdiff )
            maxdiff = diff; // enlarge maxDiff
          G[i][j] = G[i][j] + omega * (Gnew-G[i][j]); // update current cell value
        }
      }
    }
    iteration++;
  } while (maxdiff > stopdiff);

    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    printf("SOR %d x %d complete, took %ld ms to complete\n",N-2,N-2, mtime);
    printf("using %d iterations, diff is %f (allowed diff %f)\n",
           iteration,maxdiff,stopdiff);

  return 0;
}
