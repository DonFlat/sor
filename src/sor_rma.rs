use std::f64::consts::PI;
use mpi::collective::SystemOperation;
use mpi::traits::*;
use mpi::Rank;
use mpi::topology::SimpleCommunicator;
use mpi::window::{AllocatedWindow, WindowOperations};
use crate::test_utils::{append_to_csv, powers_of_two};

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

pub fn runner(problem_size: usize) {
    // ************************
    // **** Setting Up MPI ****
    // ************************
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let world_size = world.size();
    let rank = world.rank();

    let mut time_records: Vec<f64> = Vec::new();
    for _ in 0..12 {
        let time = sor(problem_size, rank, world_size, &world);
        time_records.push(time);
    }
    if rank == 0 {
        append_to_csv("rma_data.csv", problem_size, &time_records).expect("Error happened writing csv");
    }
}

pub fn sor(problem_size: usize, rank: Rank, size: Rank, world: &SimpleCommunicator) -> f64 {
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
    N += 2;

    // ************************************
    // **** Setting Up Algo parameters ****
    // ************************************
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

    // ***************************
    // **** Set up the window ****
    // ***************************
    let mut huge_window: AllocatedWindow<f64> = world.allocate_window((local_ub + 1) * n_col);
    // println!("Rank {} has {} rows", rank, local_ub+1);

    let local_ub_pred_rank = calc_pred_rank_local_ub(size, pred_rank, N);

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

    // ********************************
    // **** Start real computation ****
    // ********************************
    let t_start = mpi::time();
    let mut iteration = 0;
    loop {
        huge_window.fence();
        if pred_rank != rank {
            huge_window.put(n_col, n_col, pred_rank as usize, (local_ub_pred_rank * n_col), n_col);
        }
        if succ_rank != rank {
            huge_window.put((local_ub - 1) * n_col, n_col, succ_rank as usize, 0, n_col);
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
    return (t_end-t_start) * 1000000.0;
}

fn calc_pred_rank_local_ub(size: Rank, pred_rank: i32, mut N: usize) -> usize {
    let (mut global_lb_pred_rank, global_ub_pred_rank) = get_bounds(N - 1, size as usize, pred_rank as usize);
    if global_lb_pred_rank == 0 {
        global_lb_pred_rank = 1;
    }
    let local_ub_pred_rank = global_ub_pred_rank - (global_lb_pred_rank - 1);
    local_ub_pred_rank
}

// let mut matrix = vec![vec![0.0; n_col]; local_ub+1];
fn set_matrix(i: usize, j: usize, val: f64, matrix: &mut Vec<f64>, ncol: usize) {
    let ind = i * ncol + j;
    matrix[ind] = val;
}

fn get_matrix(i: usize, j: usize, matrix: &Vec<f64>, ncol: usize) -> f64 {
    return matrix[i * ncol + j];
}
