use std::f64::consts::PI;
use std::ffi::{c_void};
use std::os::raw::c_int;
use mpi::collective::SystemOperation;
use mpi::traits::*;
use mpi::ffi::*;
use mpi::window::{AllocatedWindow, WindowOperations};

fn even_1_odd_0(num: usize) -> usize {
    match num % 2 {
        0 => 1,
        _ => 0
    }
}

// TODO: how would you tell it is a contiguous memory block, or array of pointers?
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

pub fn sor(problem_size: usize) {
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
    let mut huge_window:AllocatedWindow<f64> = world.allocate_window(n_col * (local_ub + 1));
    println!("Rank {} has {} rows", rank, local_ub+1);

    // Initialize the boundary value
    for i in 0..=local_ub {
        for j in 0..n_col {
            if i == 0 && global_lb - 1 == 0 {
                set_matrix(i, j, 4.56, &mut huge_window.window_vector, n_col);
            } else if i == local_ub && global_ub == n_row - 1 {
                set_matrix(i, j, 9.85, &mut huge_window.window_vector, n_col);
            } else if j == 0 {
                set_matrix(i, j, 7.32, &mut huge_window.window_vector, n_col);
            } else if j == n_col - 1 {
                set_matrix(i, j, 6.88, &mut huge_window.window_vector, n_col);
            } else {
                set_matrix(i, j, 0.00, &mut huge_window.window_vector, n_col);
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

    let t_start = mpi::time();
    // Now do the real computation
    let mut iteration = 0;
    loop {
        huge_window.fence();
        if pred_rank != rank {
            unsafe {
                MPI_Put(
                    huge_window.window_vector.as_ptr().add(n_col) as *mut c_void,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    pred_rank,
                    (local_ub * n_col) as MPI_Aint,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    huge_window.window_ptr
                );
            }
        }
        if succ_rank != rank {
            unsafe {
                MPI_Put(
                    huge_window.window_vector.as_ptr().add((local_ub - 1) * n_col) as *mut c_void,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    succ_rank,
                    0,
                    n_col as c_int,
                    RSMPI_DOUBLE,
                    huge_window.window_ptr
                );
            }
        }
        huge_window.fence();
        max_diff = 0.0;
        for phase in 0..2 {
            let mut global_row_num = global_lb;
            for i in 1..local_ub {
                let start_col = 1 + (even_1_odd_0(global_row_num) ^ phase);
                for j in (start_col..n_col-1).step_by(2) {
                    // The stencil operation
                    let up = get_matrix(i - 1, j, &huge_window.window_vector, n_col);
                    let down = get_matrix(i + 1, j, &huge_window.window_vector, n_col);
                    let left = get_matrix(i, j - 1, &huge_window.window_vector, n_col);
                    let right = get_matrix(i, j + 1, &huge_window.window_vector, n_col);

                    let stencil_val = (up + down + left + right) / 4.0;
                    let mat_i_j = get_matrix(i, j, &huge_window.window_vector, n_col);
                    diff = (stencil_val as f64 - mat_i_j as f64).abs();

                    if diff > max_diff {
                        max_diff = diff;
                    }
                    set_matrix(i, j, mat_i_j + omega * (stencil_val - mat_i_j), &mut huge_window.window_vector, n_col);
                }
                global_row_num += 1;
            }
        }

        diff = max_diff;
        world.all_reduce_into(&diff, &mut max_diff, SystemOperation::max());
        iteration += 1;

        if max_diff <= stop_diff {
            break;
        }

    }
    let t_end = mpi::time();

    if rank == 0 {
        println!("SOR size: {} x {}, time: {} ms", n_row-2, n_col-2, (t_end-t_start) * 1000.0);
        println!("using {} iterations, diff is {} (allowed diff {})", iteration,max_diff,stop_diff)
    }
}

// let mut matrix = vec![vec![0.0; n_col]; local_ub+1];
fn set_matrix(i: usize, j: usize, val: f64, matrix: &mut Vec<f64>, ncol: usize) {
    let ind = i * ncol + j;
    matrix[ind] = val;
}

fn get_matrix(i: usize, j: usize, matrix: &Vec<f64>, ncol: usize) -> f64 {
    return matrix[i * ncol + j];
}
