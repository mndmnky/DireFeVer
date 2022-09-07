use fxhash::FxHashSet;
use std::error;
use crate::dfvs_instance::DFVSInstance;
use crate::reduction_rules::Rule;
use crate::digraph::Digraph;

impl DFVSInstance {

    /// Advanced branch and reduce algorithm that does the following: 
    /// 1. Exhaustive initial reduction and computation of upper and lower bounds.
    /// 1.1. The upper bound also computes an initial solution.
    /// 1.2. If the upper bound is equal to the lower bound, the initial solution is returned.
    /// 3. Priority branching:
    /// 3.1. If a cluster with more than 2 nodes is found, branch on the remaining nodes (the other
    ///   nodes are removed and added to the solution).
    /// 3.2. If a daisy with at least one petal exists, branch on the core and the petals of the
    ///   maximal daisy. 
    /// 3.3. Branch on a node with the hightest cai weight (0.3).
    /// See .advanced_branching_recursion() to see what happens in the invidual branches. 
    /// `self.changes` is reset at the start of the execution of this function.
    /// Returns the optimal solution if one was found or `None` if `self.interrupter` was set.
    ///
    /// 1.3.0 Note: With the link node rule branching on double paths has become redundent.
    ///
    /// # Arguments
    /// * `priority_rules` - A list of `Rule`s specifying which rule is applied during the
    /// execution and in which order they are applied.
    ///
    /// TODO make sure that you dont redo changes too far!!!
    pub fn advanced_branching(&mut self, priority_rules: &Vec<Rule>, skip_initial_rules: bool) -> Result<FxHashSet<usize>, Box<dyn error::Error>> {
        self._reset_changes();
        let mut best_solution;
        if !skip_initial_rules {
            // Exhaustively perform initial reductions w.o. upper bounds.
            self.exhaustive_reductions(priority_rules)
        }
        if self.graph.scc_graph_is_disconnected() {
            let mut solution = self.solution.clone();
            // create subgraph for each (keep node ids)
            let subgraphs: Vec<Digraph> = self.graph.split_into_connected_components_alt();
            for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(graph, None, None)){
                // start advanced_branching on each component
                solution.extend(instance.advanced_branching(priority_rules, true)?);
            }
            // If we use the link node rule, we have to finallize the found solution.
            solution = self.finallize_given_solution_temp(&solution);
            return Ok(solution)
        }

        // If `self.upper_bound` is `None` while `self.lower_bound` is some, we recomputed the
        // lower bound never the less.
        if self.current_best.is_none() {
            // compute (simple?) upper bound (and lower bound)
            self.compute_and_set_upper_lower(true);
        }
        if self.lower_bound.is_none() {
            self.compute_and_set_lower(true);
        }
        best_solution = self.current_best.clone().expect("is some");
        // Return solution if best upper == best lower. 
        if self.upper_bound == self.lower_bound {
            // If we use the link node rule, we have to finallize the found solution.
            best_solution = self.finallize_given_solution_temp(&best_solution);
            return Ok(best_solution)
        }
        // Start recursive branching (one recursion of .advanced_branching_recursion)
        let greedy_cluster = self.graph.greedy_max_clique();
        if greedy_cluster.len() > 2 {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            // Start with the whole cluster:
            self.add_all_to_solution(&greedy_cluster.clone()).expect("The given set of nodes exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register and sets a new (second) one.
            self.rebuild_section()?;
            self.start_new_changes();
            // Add all but one, contract the one:
            for node in &greedy_cluster {
                // Remove all but node. 
                self.add_all_to_solution(&greedy_cluster.clone().iter().copied().filter(|clu| clu !=node).collect()).expect("The given set of nodes exists");
                self.contract_node(*node).expect("This node exists");
                // Branch 
                if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                    // Return solution if best upper == best lower. 
                    // TODO: did we even change upper bound to reach this?
                    if self.upper_bound == self.lower_bound {
                        // Recovers up to very first register (pre reductions) and returns
                        self.rebuild_section()?;
                        self.rebuild_section()?;
                        return Ok(solution)
                    }
                    // Any found solution is better, than the old, due to the upper bound updates.
                    best_solution = solution;
                }
                // Recovers 2nd register and sets a new (second) one.
                self.rebuild_section()?;
                self.start_new_changes();
            }
            // Recovers up to very first register (pre reductions) and returns
            self.rebuild_section()?;
            self.rebuild_section()?;
            return Ok(best_solution)
        // Link node branch, for link nodes n which the rule couldn't opperate. 
        } else if let Some((link_node, neighs)) = self.graph.get_any_link_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_all_to_solution(&neighs.iter().copied().collect()).expect("The given set of nodes exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            self.rebuild_section()?;
            self.add_to_solution(link_node).expect("This node exists");
            self.contract_node(neighs[0]).expect("This node exists");
            self.contract_node(neighs[1]).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers up to very first register (pre reductions) and returns
            self.rebuild_section()?;
            return Ok(best_solution)
        } else if let Some((max_core, max_daisy)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_to_solution(max_core).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register.
            self.rebuild_section()?;
            self.add_all_to_solution(&max_daisy).expect("The given set of nodes exists");
            self.contract_node(max_core).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers up to very first register (pre reductions) and returns
            self.rebuild_section()?;
            return Ok(best_solution)
        } else {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            let branch_node = self.graph.get_max_weight_node(&Digraph::cai_weight, (0.3, 0f64)).expect("`self.graph` is not empty");
            self.add_to_solution(branch_node).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recovers 2nd register, contracts `node` and set a new 2nd register.
            self.rebuild_section()?;
            self.contract_node(branch_node).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    return Ok(solution)
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = solution;
            }
            // Recover to very first register.
            self.rebuild_section()?;
        }
        Ok(best_solution)
    }

    // TODO: what reductions when?
    /// Branching subroutine of the advanced branch and reduce algorithm that does the following: 
    /// 1. Exhaustively apply simple reductions
    /// 2. Check if the graph is empty
    /// 2.1. If so, check if the current solution is below the upper bound.
    /// 2.1.1. If so, return the current solution and update the upper bound. 
    /// 2.1.2. If not, return `None`.
    /// 3. Priority branching:
    /// 3.1. If a cluster with more than 2 nodes is found, branch on the remaining nodes (the other
    ///   nodes are removed and added to the solution).
    /// 3.2. If a daisy with at least one petal exists, branch on the core and the petals of the
    ///   maximal daisy. 
    /// 3.3. Branch on a node with the hightest cai weight (0.3).
    ///
    /// 1.3.0 Note: With the link node rule branching on double paths has become redundent.
    pub fn advanced_branching_recursion(&mut self, priority_rules: &Vec<Rule>) -> Result<Option<FxHashSet<usize>>, Box<dyn error::Error>> {
        let mut best_solution = None;
        // Perform some reductions.
        self.exhaustive_reductions(priority_rules);
        if self.graph.num_nodes() == 0 {
            // Check solution size
            if self.solution.len() < self.upper_bound.expect("upper bound was set") {
                // If we use the link node rule, we have to finallize the found solution.
                let sol = self.finallize_solution_temp();
                self.set_current_best(&sol);
                return Ok(Some(sol));
            } else {
                return Ok(None);
            }
        } else if self.solution.len() >= self.upper_bound.expect("upper bound was set") {
            return Ok(None);
        }
        // TODO: we could do a lower bound on all sccs and check the advanced lower bound criteria
        // we do in the next step.
        if self.graph.scc_graph_is_disconnected() {
            let mut solution = self.solution.clone();
            // create subgraph for each (keep node ids)
            let subgraphs: Vec<Digraph> = self.graph.split_into_connected_components_alt();
            for mut instance in subgraphs.into_iter().map(|graph| DFVSInstance::new(graph, None, None)){
                // start advanced_branching on each component
                solution.extend(instance.advanced_branching(priority_rules, true)?)
            }
            if solution.len() < self.upper_bound.expect("a upper bound was set") {
                // If we use the link node rule, we have to finallize the found solution.
                solution = self.finallize_given_solution_temp(&solution);
                self.set_current_best(&solution);
                return Ok(Some(solution));
            } else {
                return Ok(None);
            }
        }
        // If the current lower bound is equal or greater than our current best solution, we won't find a better
        // solution in this branch.
        let lower = self.clique_heuristic_lower(priority_rules, true);
        if lower >= self.upper_bound.expect("`upper_bound` was set") {
            return Ok(None);
        }
        // TODO (Consider): We could recompute upper bounds and check for matches with the lower
        // bound...
        let greedy_cluster = self.graph.greedy_max_clique();
        if greedy_cluster.len() > 2 {
            // In branch first register after reductions.
            self.start_new_changes();
            // Start with the whole cluster:
            self.add_all_to_solution(&greedy_cluster.clone()).expect("The given set of nodes exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recover first in branch register.
                    self.rebuild_section()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recover first in branch register and set a new one.
            self.rebuild_section()?;
            self.start_new_changes();
            // Add all but one, contract the one:
            for node in &greedy_cluster {
                // Remove all but node. 
                self.add_all_to_solution(&greedy_cluster.clone().iter().copied().filter(|clu| clu != node).collect()).expect("The given set of nodes exists");
                self.contract_node(*node).expect("This node exists");
                // Branch 
                if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                        // Return solution if best upper == best lower. 
                        if self.upper_bound == self.lower_bound {
                            // Recover first in branch register.
                            self.rebuild_section()?;
                            return Ok(Some(solution))
                        }
                    // Any found solution is better, than the old, due to the upper bound updates.
                    best_solution = Some(solution);
                }
                // Recover first in branch register and set new one.
                self.rebuild_section()?;
                self.start_new_changes();
            }
            // Recover first in branch register.
            self.rebuild_section()?;
            return Ok(best_solution)
        // Link node branch, for link nodes n which the rule couldn't opperate. 
        } else if let Some((link_node, neighs)) = self.graph.get_any_link_node_and_neighbors() {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            self.add_all_to_solution(&neighs.iter().copied().collect()).expect("The given set of nodes exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            self.rebuild_section()?;
            self.add_to_solution(link_node).expect("This node exists");
            self.contract_node(neighs[0]).expect("This node exists");
            self.contract_node(neighs[1]).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recovers up to very first register (pre reductions) and returns
            return Ok(best_solution)
        } else if let Some((max_core, max_daisy)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            self.start_new_changes();
            self.add_to_solution(max_core).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // End of for loop to pop last register
                    self.rebuild_section()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            self.rebuild_section()?;
            self.add_all_to_solution(&max_daisy).expect("The given set of nodes exists");
            self.contract_node(max_core).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            return Ok(best_solution)
        } else {
            // 2nd register (after initial reductions):
            self.start_new_changes();
            let branch_node = self.graph.get_max_weight_node(&Digraph::cai_weight, (0.3, 0f64)).expect("`self.graph` is not empty");
            self.add_to_solution(branch_node).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    self.rebuild_section()?;
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recovers 2nd register, contracts `node` and set a new 2nd register.
            self.rebuild_section()?;
            self.contract_node(branch_node).expect("This node exists");
            if let Some(solution) = self.advanced_branching_recursion(priority_rules)? {
                // Return solution if best upper == best lower. 
                if self.upper_bound == self.lower_bound {
                    // Recovers up to very first register (pre reductions) and returns
                    return Ok(Some(solution))
                }
                // Any found solution is better, than the old, due to the upper bound updates.
                best_solution = Some(solution);
            }
            // Recover to very first register.
        }
        Ok(best_solution)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::digraph::Digraph;
    use crate::dfvs_instance::DFVSInstance;

    #[test]
    fn advanced_branch_test() {
        let gr = Cursor::new("48 138 0\n2 3 4 5 12\n1 3 4 5 6 7 8\n1 2 4 5 14 15 16\n\
                             1 2 3 5 9 10 11\n1 2 3 4 9 10 11\n2 7\n2 8\n2 6 13\n\
                             4 5\n4 5\n4 5\n1 13\n6 12\n3 15\n3 16\n3 14\n18 19 20 21\n\
                             17 19\n17 20\n17 21\n17 18\n23 24 25\n22 24 26 28\n22 23 29 30 31\n\
                             22 27\n23 27\n23 28\n25 26\n24 30\n24 31\n24 29\n33 34 35 36\n\
                             32 37 38 39\n32 40 41 42\n32 43 44 45\n32 46 47 48\n\
                             33 38\n33 39\n33 37\n34 41\n34 42\n34 40\n35 44\n35 45\n\
                             35 43\n36 47\n36 48\n36 46\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let graph = g;
        let mut branch_instance = DFVSInstance::new(graph, None, None);
        let solution = branch_instance.advanced_branching(&vec![Rule::SimpleRules],false);
        assert!(solution.is_ok());
        assert_eq!(solution.unwrap().len(), 22); 
    }
}
