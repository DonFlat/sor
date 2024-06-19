#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"

void write_to_csv(double window_size, int node_num, double* latency) {
    // Open a file for appending
    FILE *fpt = fopen("sor_lock_c.csv", "a");
    if (fpt == NULL) {
        printf("Error opening the file.\n");
        return;
    }
    // Write the window_size as the first column
    fprintf(fpt, "%d", (int)window_size);
    fprintf(fpt, ",%d", node_num);

    // Write the elements of the latency array as the rest of the columns
    for (int i = 0; i < 12; i++) {
        fprintf(fpt, ",%f", latency[i]);
    }
    // End the line for CSV row
    fprintf(fpt, "\n");

    // Close the file
    fclose(fpt);

    printf("Data appended to sor_c_rma.csv\n");
}

int even(int i){
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
static void get_bounds(int *lb, int *ub, int n, int size, int rank)
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

double run_sor(int size, int rank, int argc, char *argv[]) {
  int N;                            /* problem size */
  int ncol,nrow;                    /* number of rows and columns */
  double TOLERANCE = 0.0002;        /* termination criterion */
  double Gnew,r,omega,              /* some float values */
         stopdiff,maxdiff,diff;     /* differences btw grid points in iters */
  double **G;                       /* the grid */
  int i,j,phase,iteration;          /* counters */
  double t_start,t_end;             /* time values */
  
  int lb,ub;                        /* lower and upper bound of grid stripe */
  MPI_Status status;                /* dummy for MPI_Recv */


  /* ranks of predecessor and successor for row exchanges */
  int pred = (rank == 0) ? MPI_PROC_NULL : rank - 1;
  int succ = (rank == size-1) ? MPI_PROC_NULL : rank + 1;

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

  /* initialize the window */
    MPI_Win window_lb_1;
    MPI_Win_create(
        G[lb - 1],
        ncol * sizeof(double),
        sizeof(double),
        MPI_INFO_NULL,
        MPI_COMM_WORLD,
        &window_lb_1
    );
    
    MPI_Win window_ub;
    MPI_Win_create(
        G[ub],
        ncol * sizeof(double),
        sizeof(double),
        MPI_INFO_NULL,
        MPI_COMM_WORLD,
        &window_ub
    );

    /* now do the "real" computation */
    t_start = MPI_Wtime();
    iteration = 0;
    do {
        if (rank != 0) {
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, pred, 0, window_ub);
            MPI_Put(G[lb], ncol, MPI_DOUBLE, pred, 0, ncol, MPI_DOUBLE, window_ub);
            MPI_Win_unlock(pred, window_ub);
        }
        if (rank != size - 1) {
          MPI_Win_lock(MPI_LOCK_EXCLUSIVE, succ, 0, window_lb_1);
          MPI_Put(G[ub - 1], ncol, MPI_DOUBLE, succ, 0, ncol, MPI_DOUBLE, window_lb_1);
          MPI_Win_unlock(succ, window_lb_1);
        }
      MPI_Barrier(MPI_COMM_WORLD);

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

    MPI_Win_free(&window_lb_1);
    MPI_Win_free(&window_ub);

    return (t_end-t_start) * 1000000;
}

int main (int argc, char *argv[]){

  int size,rank;         /* process ranks */
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double time_costs[12];
  double time_cost = 0;
  for (int i = 0; i < 12; i++) {
    time_cost = run_sor(size, rank, argc, argv);
    time_costs[i] = time_cost;
  }
  if (rank == 0) {
    write_to_csv(atoi(argv[1]), atoi(argv[2]), time_costs);
  }
  MPI_Finalize();
  return 0;
}