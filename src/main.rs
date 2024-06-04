use std::env;

mod sor_sendrecv;
mod sor_rma_huge_window;
mod sor_seq;
mod test_utils;
mod sor_rma_raw;
mod sor_rma_split_matrix;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mpi_type: &str = &args[1];
    let problem_size: usize = args[2].parse().expect("Failed to parse args[2] as usize");
    let node_num: usize = args[3].parse().expect("Failed to parse args[3] as usize");
    match mpi_type {
        "big" => sor_rma_huge_window::runner(problem_size, node_num),
        "split" => sor_rma_split_matrix::runner(problem_size, node_num),
        "raw" => sor_rma_raw::runner(problem_size, node_num),
        "norm" => sor_sendrecv::runner(problem_size, node_num),
        _ => println!("Invalid arguments"),
    }
}
