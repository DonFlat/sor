use std::error::Error;
use std::fs::OpenOptions;
use csv::WriterBuilder;

// Leiden has max 18 nodes, 2^5 = 32. To prevent case that insufficient row for a node happen
pub fn powers_of_two(n: u32) -> Vec<u32> {
    (5..n).map(|i| 2_u32.pow(i)).collect()
}

pub fn append_to_csv(file_path: &str, vector_size: usize, repetitions: &Vec<f64>) -> Result<(), Box<dyn Error>> {
    let file = OpenOptions::new()
        .write(true)
        .append(true)
        .create(true)
        .open(file_path)?;

    let mut writer = WriterBuilder::new()
        .has_headers(false)
        .from_writer(file);

    let mut row = vec![vector_size.to_string()];
    row.extend(repetitions.iter().map(|&val| val.to_string()));

    writer.write_record(&row)?;
    writer.flush()?;

    Ok(())
}