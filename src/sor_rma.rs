use std::f64::consts::PI;
use std::ffi::{c_double, c_void};
use std::mem::size_of;
use std::os::raw::c_int;
use std::{env, ptr};
use mpi::collective::SystemOperation;
use mpi::ffi::*;
use mpi::Rank;
use mpi::traits::*;

fn even_1_odd_0(num: usize) -> usize {
    match num % 2 {
        0 => 1,
        _ => 0
    }
}

fn get_bounds(n: usize, size: usize, rank: usize) -> (usize, usize) {
    let mut nlarge = n % size; // 1000 % 4 = 0
    let mut nsmall = size - nlarge; // 4 - 0 = 4

    let mut size_small = n / size; // 1000 / 4 = 25
    let  size_large = size_small + 1; // 25 + 1 = 26

    let mut lower_bound;
    let mut upper_bound;

    if rank < nlarge { // 2 < 0 ?
        lower_bound = rank * size_large;
        upper_bound = lower_bound + size_large;
    } else {
        // 0 * 26 + (2 - 0) * 4 = 8
        lower_bound = nlarge * size_large + (rank - nlarge) * size_small;
        // 8 + 25 = 33
        upper_bound = lower_bound + size_small;
    }
    (lower_bound, upper_bound)
}

pub fn sor(problem_size: usize, iteration_limit: usize) {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let size = world.size();
    let rank = world.rank();

    let pred_rank = if rank == 0 { 0 } else { rank - 1 };
    let succ_rank = if rank == size - 1 { rank } else { rank + 1 };

    let mut N: usize  = problem_size;
    if N == 0 {
        N = 1000;
    }
    // Give each process at least one row
    if N < size as usize {
        N = size as usize;
    }

    if rank == 0 {
        println!("Running SOR on nodes: {}, rows: {}", size, N);
    }
    N += 2;

    let n_col = N;
    let n_row = N;
    let r = 0.5 * ((PI / n_col as f64).cos() + (PI / n_row as f64).cos());
    let mut  omega =  2.0 / (1.0 + (1.0 - r * r).sqrt());
    let stop_diff = 0.0002 / (2.0 - omega);
    let mut max_diff;
    let mut diff;
    omega *= 0.8;

    // get my stripe bounds and malloc the grid accordingly
    let (mut global_lb, global_ub) = get_bounds(N - 1, size as usize, rank as usize);
    // row 0 is static
    if global_lb == 0 {
        global_lb = 1;
    }

    let local_ub = global_ub - (global_lb - 1);
    let mut matrix = vec![vec![0.0; n_col]; local_ub+1];

    // Initialize the boundary value
    for i in 0..=local_ub {
        for j in 0..n_col {
            if i == 0 && global_lb - 1 == 0 {
                matrix[i][j] = 4.56;
            } else if i == local_ub && global_ub == n_row - 1 {
                matrix[i][j] = 9.85;
            } else if j == 0 {
                matrix[i][j] = 7.32;
            } else if j == n_col - 1 {
                matrix[i][j] = 6.88;
            } else {
                matrix[i][j] = 0.0;
            }
        }
    }

    // cur rank --- matrix[1] ---> pred rank matrix[local_ub]
    // cur rank --- matrix[local_ub-1] ---> succ rank matrix[0]

    // cur rank matrix[0] <-- matrix[local_ub-1] -- pred rank
    // cur rank matrix[local_ub] <-- matrix[1] -- succ rank

    // open up local_ub and 0 as windows,
    // let cur put matrix[1] to pred's local ub window
    // let cur put matrix[local_ub-1] to succ's matrix[0]
    let mut window_local_ub = ptr::null_mut();
    unsafe {
        MPI_Win_create(
            matrix[local_ub].as_mut_ptr() as *mut c_void,
            (n_col * size_of::<c_double>()) as MPI_Aint,
            size_of::<c_double>() as c_int,
            RSMPI_INFO_NULL,
            RSMPI_COMM_WORLD,
            &mut window_local_ub
        );
    }

    let mut window_row_0 = ptr::null_mut();
    unsafe {
        MPI_Win_create(
            matrix[0].as_mut_ptr() as *mut c_void,
            (n_col * size_of::<c_double>()) as MPI_Aint,
            size_of::<c_double>() as c_int,
            RSMPI_INFO_NULL,
            RSMPI_COMM_WORLD,
            &mut window_row_0
        );
    }

    let t_start = mpi::time();
    // Now do the real computation
    let mut iteration = 0;
    loop {
        unsafe {
            MPI_Win_fence(0, window_local_ub);
            MPI_Win_fence(0, window_row_0);
        }
        if pred_rank != rank {
            unsafe {
                MPI_Put(
                    matrix[1].as_mut_ptr() as *mut c_void,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    pred_rank,
                    0,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    window_local_ub
                );
            }
        }
        if succ_rank != rank {
            unsafe {
                MPI_Put(
                    matrix[local_ub-1].as_mut_ptr() as *mut c_void,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    succ_rank,
                    0,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    window_row_0
                );
            }
        }
        unsafe {
            MPI_Win_fence(0, window_local_ub);
            MPI_Win_fence(0, window_row_0);
        }
        max_diff = 0.0;
        for phase in 0..2 {
            let mut global_row_num = global_lb;
            for i in 1..local_ub {
                let start_col = 1 + (even_1_odd_0(global_row_num) ^ phase);
                for j in (start_col..n_col-1).step_by(2) {
                    // The stencil operation
                    let up = &matrix[i - 1][j];
                    let down = &matrix[i + 1][j];
                    let left = &matrix[i][j - 1];
                    let right = &matrix[i][j + 1];

                    let stencil_val = (up + down + left + right) / 4.0;
                    diff = (stencil_val as f64 - matrix[i][j] as f64).abs();

                    if diff > max_diff {
                        max_diff = diff;
                    }
                    matrix[i][j] = matrix[i][j] + omega * (stencil_val - matrix[i][j]);
                }
                global_row_num += 1;
            }
        }

        diff = max_diff;
        // TODO: can remove it to check pure rma, fixed iterations
        world.all_reduce_into(&diff, &mut max_diff, SystemOperation::max());
        iteration += 1;

        if iteration == iteration_limit {
            break;
        }

    }
    let t_end = mpi::time();

    if rank == 0 {
        println!("SOR size: {} x {}, time: {} ms", n_row-2, n_col-2, (t_end-t_start) * 1000.0);
        println!("using {} iterations, diff is {} (allowed diff {})", iteration,max_diff,stop_diff)
    }
}

fn print_matrix(matrix: &Vec<Vec<f64>>, n_row: usize, n_col: usize, rank: Rank) {
    println!("Rank {}'s matrix -----------", rank);
    for i in 0..n_row {
        print!("r[{}]: ", i);
        for j in 0..n_col {
            print!(" {} ", matrix[i][j]);
        }
        println!();
    }
}

// TODO: revisit RMA paper, why small message is quicker?
// TODO: try ping-pong test
// TODO: SOR with different numbers of iterations, problem size. don't need to get correct result
// Some found for experiment with: 8 nodes, 10 iters,
// matrix size: 50.000, RMA: 8460.146122999999 ms, sendrev too long to get result
// Even C sor -O3 on 50K took too long to get result, but 30k took 3.166942 s, RMA is 3562.165272 ms
// matrix size: 40.000, RMA: 5720.041005 ms, sendrev also too long to get result, C sor is same
//
// TODO: seq C vs seq Rust

// TODO: library, refactor test, revisit RMA paper, know synchronization
// TODO: see CPU cycle
// start from simple case regardless of displacement unit or so.