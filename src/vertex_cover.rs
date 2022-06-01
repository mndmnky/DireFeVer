use duck_and_cover::vc_instance::VCInstance;
use duck_and_cover::graph::DyUGraph;
use crate::{digraph::Digraph, cust_errors::ImportError};
use std::io::prelude::*;
use fxhash::FxHashSet;

impl From<Digraph> for VCInstance {
    fn from(graph: Digraph) -> Self {
        let dyugraph = DyUGraph::from_edge_iter(graph.strong_edges(), graph.num_reserved_nodes());
        VCInstance::new(dyugraph)
    }
}

/// Read a vertex cover solution and returns it.
pub fn read_vc_solution<R: BufRead>(gr: R) -> Result<FxHashSet<usize>, ImportError> {
        let (lines, _comments): (Vec<_>, Vec<_>) = gr.lines()
            .partition(|l| {
                if let Ok(line) = l {
                    // separate empty lines and comment lines
                    !line.starts_with("c ")
                } else {
                    true
                }
            });
        let mut solution = FxHashSet::default();
        let mut lines = lines.into_iter();
        let _s = {
            let line = lines.next().ok_or(ImportError::InputMalformedError)??;
            let mut split = line.split(' ');
            if let Some("s") = split.next() {} else { return Err(ImportError::InputMalformedError); }
            if let Some("vc") = split.next() {} else { return Err(ImportError::InputMalformedError); }
            let _num_nodes: usize = split.next().ok_or(ImportError::InputMalformedError)?.parse()?;
            let s: usize = split.next().ok_or(ImportError::InputMalformedError)?.parse()?;
            if split.next().is_some() { return Err(ImportError::InputMalformedError); }
            s
        };
        for line in lines {
            let sol_vert: usize = line?.parse()?;
            solution.insert(sol_vert);
        }
        Ok(solution)
}
