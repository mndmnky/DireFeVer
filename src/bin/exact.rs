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
    let priority = vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
    dfvsi.exhaustive_reductions(&priority);
    let mut resu;
    if dfvsi.graph.scc_graph_is_disconnected() {
        let mut solution = dfvsi.solution.clone();
        // create subgraph for each (keep node ids)
        let subgraphs: Vec<Digraph> = dfvsi.graph.split_into_connected_components_alt();
        for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(graph, None, None)){
            // check size and other important attributes to see if it makes sense to switch to vc.
            if instance.graph.is_vcable() {
                let i = 0;
                match instance.via_vertex_cover() {
                    VCQ::Exact(sol) => {
                        // If vc is a solution, add to solution, 
                        solution.extend(&sol);
                        continue
                    },
                    VCQ::Bounds(lower, opt_upper) => {
                        // Else, set lower bound
                        instance.update_lower_bound(lower);
                        instance.set_current_best(&opt_upper);
                    }
                }
            }
            // Do advanced_branching and add solution 
            let sol = instance.advanced_branching(&priority, true)?;
            solution.extend(&sol);
        }
        // If we use the link node rule, we have to finallize the found solution.
        resu = dfvsi.finallize_given_solution_temp(&solution);
    } else {
        if dfvsi.graph.is_vcable() {
            match dfvsi.via_vertex_cover() {
                VCQ::Exact(sol) => {
                    // If vc is a solution, add to solution, 
                    // finallize
                    resu = dfvsi.solution.clone();
                    resu.extend(&sol);
                    resu = dfvsi.finallize_given_solution_temp(&resu);
                },
                VCQ::Bounds(lower, opt_upper) => {
                    // Else, set lower bound
                    dfvsi.update_lower_bound(lower);
                    dfvsi.set_current_best(&opt_upper);
                    resu = dfvsi.advanced_branching(&priority, true)?;
                }
            }
        } else {
            // Do advanced_branching and add solution 
            resu = dfvsi.advanced_branching(&priority, true)?;
        }
    }
    DFVSInstance::write_solution(&resu, &mut stdout)?;
    Ok(())
}
