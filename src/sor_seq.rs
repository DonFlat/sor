use std::{env, time};
use std::f64::consts::PI;
use std::time::Instant;

fn stencil(matrix: &Vec<Vec<f64>>, row: usize, col: usize) -> f64 {
    (matrix[row - 1][col] + matrix[row + 1][col] + matrix[row][col - 1] + matrix[row][col + 1]) / 4.0
}

fn even_1_odd_0(num: usize) -> usize {
    match num % 2 {
        0 => 1,
        _ => 0
    }
}
pub fn sor(problem_size: usize, iteration_limit: usize) {
    // Problem size, default 1K (the matrix size)
    let mut N: usize  = problem_size;
    if N == 0 {
        N = 1000;
    }

    println!("Running SOR with {} rows", N);

    // Adding two borderlines
    N += 2;

    // Set up matrix info
    let n_col = N;
    let n_row = N;

    // Initialize parameters
    let r = 0.5 * ((PI / n_col as f64).cos() + (PI / n_row as f64).cos());
    let mut  omega =  2.0 / (1.0 + (1.0 - r * r).sqrt());
    let stop_diff = 0.0002 / (2.0 - omega);
    let mut max_diff;
    let mut diff;
    omega *= 0.8;

    // Initialize the matrix
    let mut mat = vec![vec![0.0; n_col]; n_row];
    let mut m = [[0.0; 1000]; 1000];
    let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(n_row);
    for _ in 0..n_row {
        matrix.push(vec![0.0; n_col]);
    }
    for i in 0..n_row {
        for j in 0..n_col {
            if i == 0 {
                matrix[i][j] = 4.56;
            } else if i == n_row - 1 {
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

    // Now do the real computation
    let start = Instant::now();
    let mut iteration: i32 = 0;
    loop {
        max_diff = 0.0;
        for phase in 0..2 {
            for i in 1..(N-1) {
                let start_col = 1 + (even_1_odd_0(i) ^ phase);
                for j in (start_col..(N-1)).step_by(2) {
                    let stencil_val = stencil(&matrix, i, j);
                    diff = (stencil_val - matrix[i][j]).abs();
                    if diff > max_diff {
                        max_diff = diff;
                    }
                    matrix[i][j] = matrix[i][j] + omega * (stencil_val - matrix[i][j]);
                }
            }
        }
        iteration += 1;
        if max_diff <= stop_diff {
            break;
        }
    }
    let elapsed = start.elapsed();

    // Print results
    println!("SOR size: {} x {}, time: {} ms", n_row-2, n_col-2, elapsed.as_millis());
    println!("Using {} iterations, diff is {} (allowed diff {})", iteration, max_diff, stop_diff);
}