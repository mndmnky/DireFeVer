//!
//! This binary is only meant for experiments.
//! Writes some basic stats of the input graphs to a csv file. Those stats are: 
//! * Number of nodes 
//! * Number of edges 
//! * Number of PIE-edges
//! * An simple upper bound
//! * The time it took to compute the upper bound 
//! * A simple lower bound
//! * The time it took to compute the lower bound 

use std::error;
use std::time::Instant;
use std::sync::mpsc::SendError;
use clap::{Arg, Command};
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufReader, Write};
use std::fmt::Display;

use dfvs_solver::digraph::Digraph;
use dfvs_solver::dfvs_instance::DFVSInstance;

#[derive(Debug)]
enum ThreadErr {
    SendError(SendError<i8>),
    IoError(std::io::Error),
}

impl Display for ThreadErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThreadErr::SendError(send_error) => 
                write!(f, "{}", send_error),
            ThreadErr::IoError(io_error) => 
                write!(f, "{}", io_error),
        }
    }
}

impl std::error::Error for ThreadErr {}

impl From<SendError<i8>> for ThreadErr {
    fn from(err: SendError<i8>) -> Self {
        ThreadErr::SendError(err)
    }
}

impl From<std::io::Error> for ThreadErr {
    fn from(err: std::io::Error) -> Self {
        ThreadErr::IoError(err)
    }
}

pub fn main() -> Result<(), Box<dyn error::Error>> {
    // CLI stuff
    let m = Command::new("init")
        .arg(Arg::new("files")
             .takes_value(true)
             .multiple_values(true)
             .short('f'))
        .arg(Arg::new("dest")
             .required(true)
             .takes_value(true)
             .short('d'))
        .get_matches();
    // Get, as input all public instances. 
    let files: Vec<PathBuf> = m.values_of("files").unwrap().map(|p| PathBuf::from(p)).collect();
    let dest: &str = m.value_of("dest").unwrap();

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/stats.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, n, m, pie_edges, upper, t_upper, lower, t_lower")?;

    // Read graphs
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let mut dfvsi = DFVSInstance::new(graph, None, None);
        let name = file.file_stem().expect("Not a file.");
        let pie_edges_count = dfvsi.graph.strong_edges().count();
        let start_time = Instant::now();
        dfvsi.compute_and_set_fast_upper(false);
        let time_upper = start_time.elapsed().as_millis();
        let start_time = Instant::now();
        dfvsi.compute_and_set_lower(false);
        let time_lower = start_time.elapsed().as_millis();
        let upper_init = dfvsi.upper_bound.expect("was set");
        let lower_init = dfvsi.lower_bound.expect("was set");
        writeln!(out_files[0], "{:?}, {}, {}, {}, {}, {}, {}, {}", name, graph.num_nodes(), 
                 graph.num_edges(), pie_edges_count, upper_init, time_upper, lower_init, time_lower)?;
    }
    Ok(())
}

