use std::env;

mod sor_sendrecv;
mod sor_rma;
mod sor_seq;
mod test_utils;
mod sor_rma_raw;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mpi_type: &str = &args[1];
    let size: usize = args[2].parse().expect("Failed to parse args[2] as u32");
    match mpi_type {
        "rma" => sor_rma::runner(size),
        "raw" => sor_rma_raw::runner(size),
        "norm" => sor_sendrecv::runner(size),
        _ => println!("Invalid argument, run either ping pong | sor_source_data, rma | norm"),
    }
}
