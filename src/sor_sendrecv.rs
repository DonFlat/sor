use std::env;
use std::f64::consts::PI;
use mpi::collective::SystemOperation;
use mpi::Rank;
use mpi::topology::SimpleCommunicator;
use mpi::traits::*;
use crate::test_utils::powers_of_two;

fn even_1_odd_0(num: usize) -> usize {
    match num % 2 {
        0 => 1,
        _ => 0
    }
}

// Suppose n = 1000, size = 4, rank = 2
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
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();

    for n in powers_of_two(problem_size as u32) {
        sor(n as usize, rank, &world);
    }
}

pub fn sor(problem_size: usize, rank: Rank, world: &SimpleCommunicator) {
    let size = world.size();

    let pred_rank = if rank == 0 { rank } else { rank - 1 };
    let succ_rank = if rank == size - 1 { rank } else { rank + 1 };

    // println!("I'm rank {}, my pred rank: {}, succ rank: {}", rank, pred_rank, succ_rank);
    // initialize the basic variables
    let mut N = problem_size;
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
    // println!("Stop diff: {}", stop_diff);
    let mut max_diff;
    let mut diff;
    omega *= 0.8;

    // get my stripe bounds and malloc the grid accordingly
    let (mut global_lb, mut global_ub) = get_bounds(N - 1, size as usize, rank as usize);
    // row 0 is static
    if global_lb == 0 {
        global_lb = 1;
    }
    // Initialize the matrix at local rank, full size, 0 filled
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

    let t_start = mpi::time();
    // Now do the real computation
    let mut iteration = 0;
    loop {
        if pred_rank != rank {
            world.process_at_rank(pred_rank).send_with_tag(&matrix[1], 42);
            world.process_at_rank(pred_rank).receive_into_with_tag(&mut matrix[0], 42);
        }
        if succ_rank != rank {
            world.process_at_rank(succ_rank).send_with_tag(&matrix[local_ub -1], 42);
            world.process_at_rank(succ_rank).receive_into_with_tag(&mut matrix[local_ub], 42);
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
// try from 10x10, check matrix if identical
// iterations should be same
// only use boundary rows
// check out if optimization options