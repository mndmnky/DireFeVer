//! bin for quick tests
//!
//! Test how much semi node disjunct cycles are found with different parameters

use std::error;
use std::io;

use dfvs_solver::{digraph::Digraph, dfvs_instance::DFVSInstance, reduction_rules::Rule, heuristics::VCQ};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let stdin = io::stdin();
    let stdin = stdin.lock();
    let stdout = io::stdout();
    let mut stdout = stdout.lock();
    let graph = Digraph::read_graph(stdin)?;
    let mut dfvsi = DFVSInstance::new(graph, None, None);
    let priority = vec![Rule::SimpleRules];
    dfvsi.exhaustive_reductions(&priority);
    for c in 1..7 {
        let cycles = dfvsi.graph.find_semi_disjoint_cycles(c);
        eprintln!("c = {}",c);
        eprintln!("cycles = {}", cycles.len());
    }
    Ok(())
}
