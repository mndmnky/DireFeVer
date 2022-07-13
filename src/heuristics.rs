use crate::dfvs_instance::DFVSInstance;
use duck_and_cover::{vc_instance::VCInstance, kernelization::RECOMMENDED};
use fxhash::FxHashSet;
use crate::digraph::Digraph;
use crate::reduction_rules::Rule;

pub enum VCQ {
    Exact(FxHashSet<usize>),
    Bounds(usize,FxHashSet<usize>),
}

impl DFVSInstance {

    /// Converts `self.graph.graph` into a vertex cover instance, where only strong edges are
    /// considered, and solves it.
    /// Checks if the found solution is also a solution for the `DFVSInstance`
    ///
    /// `vc_instance` can not be interrupted.
    pub fn via_vertex_cover(&self) -> VCQ {
        let mut vc_instance: VCInstance = self.graph.clone().into();
        let resu = vc_instance.branch_and_reduce(RECOMMENDED).expect("Should not throw a `ProcessingError` on a fresh instance");
        return if self.validate_left_over(&resu) {
            VCQ::Exact(resu)
        } else {
            let lower = resu.len();
            // also get remaining upperbound.
            let mut clone = self.clone();
            clone.add_all_to_solution(&resu).expect("All those nodes exist.");
            let upper = clone.get_good_upper(false);
            VCQ::Bounds(lower, upper)
        }
    }

    /// Exhaustively applies local search on a heuristic (complete) solution 
    /// `some_solution` (but only actually considers the delta set of the solution that was added
    /// by the heuristic).
    ///
    /// Returns a better solution if one was found.
    pub fn exhaustive_local_search(&self, some_solution: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        // Find delta solution 
        let delta_solution: FxHashSet<usize> = some_solution.difference(&self.solution).copied().collect();
        let mut new_solution = delta_solution.clone();
        let mut changed = false;
        for node in delta_solution {
            // Can there be nodes that aren't in graph
            if !self.graph.has_node(node) {
                continue;
            }
            let mut graph_clone = self.graph.clone();
            let mut reduced_set = new_solution.clone();
            reduced_set.remove(&node);
            graph_clone.remove_nodes(reduced_set.clone());
            if !graph_clone.has_cycle() {
                changed = true;
                new_solution = reduced_set;
            }
        }
        // Rebuild solution and return
        if changed {
            let new_solution: FxHashSet<_> = new_solution.union(&self.solution).copied().collect();
            return Some(new_solution)
        }
        return None
    }

    /// Applies a local search over some heuristic (delta) solution with the goal to 
    /// remove a redundent node. 
    /// Returns some `new_solution` with |`new_solution`| = |`delta_solution`| - 1, or `None`
    /// otherwise.
    pub fn local_search(&self, delta_solution: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        for node in delta_solution {
            let mut graph_clone = self.graph.clone();
            let mut reduced_set = delta_solution.clone();
            reduced_set.remove(&node);
            graph_clone.remove_nodes(reduced_set.iter().cloned());
            if !graph_clone.has_cycle() {
                return Some(reduced_set)
            }
        }
        None
    }

    /// Performs a top down weight heuristic, where repeatedly the node with the highest weight, computed by
    /// `fn` given `weights`, is removed and added to the solution and some reduction rules are applied.
    /// At each step each rule in `rule_priority` is executed exhaustively as described in
    /// `self.exhaustive_reductions()`.
    /// Returns the solution found by the heuristic.
    ///
    /// # Panics
    /// Panics if the first rule of `rule_priority` is not `Rule::SimpleRules`.
    pub fn top_down_weight_heuristic(&self, fun: &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64), rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> FxHashSet<usize> {
        let mut clone_instance = self.clone();
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            if let Some(hw_node) = clone_instance.graph.get_max_weight_node(fun, weights) {
                clone_instance.solution.insert(hw_node);
                clone_instance.graph.remove_node(hw_node);
            } else {
                clone_instance.finallize_solution();
                return clone_instance.solution
            }
        }
    }

    /// Performs a top/bottom weight heuristic, where repeatedly the node with the lowest weight, computed by
    /// `fn` given `weights`, is contracted, some reduction rules are applied and then the node
    /// with the highest weight is removed, added to the solution and the remaining instance is
    /// reduced again.
    /// At each step each rule in `rule_priority` is executed exhaustively as described in
    /// `self.exhaustive_reductions()`.
    /// Returns the solution found by the heuristic.
    ///
    /// # Panics
    /// Panics if the first rule of `rule_priority` is not `Rule::SimpleRules`.
    pub fn top_bottom_switch_weight_heuristic(&self, fun: &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64), rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> FxHashSet<usize> {
        let mut clone_instance = self.clone();
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            if let Some(hw_node) = clone_instance.graph.get_min_weight_node(fun, weights) {
                clone_instance.graph.contract_node(hw_node);
            } else {
                clone_instance.finallize_solution();
                return clone_instance.solution
            }
            // perform reductions:
            clone_instance.exhaustive_reductions(rule_priority);
            if let Some(hw_node) = clone_instance.graph.get_max_weight_node(fun, weights) {
                clone_instance.solution.insert(hw_node);
                clone_instance.graph.remove_node(hw_node);
            } else {
                clone_instance.finallize_solution();
                return clone_instance.solution
            }
        }
    }

    /// Performs a bottom up weight heuristic, where repeatedly the node with the lowest weight, computed by
    /// `fn` given `weights`, is contracted and some reduction rules are applied.
    ///
    /// Returns the solution found by the heuristic.
    ///
    /// # Panics
    /// Panics if the first rule of `rule_priority` is not `Rule::SimpleRules`.
    pub fn bottom_up_weight_heuristic(&self, fun: &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64), rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> FxHashSet<usize> {
        let mut clone_instance = self.clone();
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            if let Some(hw_node) = clone_instance.graph.get_min_weight_node(fun, weights) {
                clone_instance.graph.contract_node(hw_node);
            } else {
                clone_instance.finallize_solution();
                return clone_instance.solution
            }
        }
    }

    /// Heuristic that does the following: 
    /// 1. Apply `rule_priority` exhaustively, or skips it if `skip_initial_rules` is set to
    ///    `true`.
    /// 2. Greedily find an inclusion maximal clique `cliq`
    /// 3. Removes all nodes in `cliq` from the graph
    /// 4. Adds `cliq`.len() to the upper bound and `cliq`.len() - 1 to the lower bound
    /// 5. Repeat from 1. until no more clique is found
    /// 6. The solution (that is found by the reduction rules) size is added to both bounds.
    /// 7. Add the upper bound from `.highest_degree_heuristic()` to the upper bound.
    /// 8. Add the lower bound from `.graph.small_disjunct_cycles_heuristic()` to the lower bound 
    /// TODO: Could use `local_search()`
    ///
    /// Returns a lower- and a upper bound and the solution of the upper bound.
    ///
    /// # Panics
    /// Panics if the first rule of `rule_priority` is not `Rule::SimpleRules`.
    pub fn clique_heuristic(&self, rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> (usize, usize, FxHashSet<usize>) {
        let mut clone_instance = self.clone();
        let mut lower_bound = 0;
        let mut upper_bound = 0;
        let mut upper_sol: FxHashSet<usize> = FxHashSet::default();
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert!(!rule_priority.contains(&Rule::AdvancedPetal));
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            // Greedily find clique `cliq`
            let cliq = clone_instance.graph.greedy_max_clique();
            if cliq.len() > 1 {
                // TODO: remove cliques with higher degree first (especially when `cliq.len() ==
                // 2`.
                // Remove `cliq` from graph 
                clone_instance.graph.remove_nodes(cliq.clone());
                lower_bound += cliq.len() - 1;
                upper_bound += cliq.len();
                upper_sol.extend(cliq);
            } else {
                // until no more clique is found.
                break
            }
        }
        // add solution size to both lower and upper.
        lower_bound += clone_instance.solution.len();
        upper_bound += clone_instance.solution.len();
        clone_instance.finallize_solution();
        upper_sol.extend(clone_instance.solution);
        // find bounds for left overs
        let instance = DFVSInstance::new(clone_instance.graph.clone(), None, None);
        // TODO: adapt to real best upper heuristic.
        let add_upper_sol = instance.top_down_weight_heuristic(&Digraph::cai_weight, (0.2,0.0), rule_priority, true) ;
        upper_bound += add_upper_sol.len();
        upper_sol.extend(add_upper_sol);
        // TODO: could be improved.
        let add_lower_sol = instance.disjunct_cycle_heuristic(rule_priority, true);
        lower_bound += add_lower_sol;
        let upper_sol = instance.finallize_given_solution_temp(&upper_sol);
        return (lower_bound, upper_bound, upper_sol)
    }

    /// Heuristic that does the following: 
    /// 1. Apply `rule_priority` exhaustively, or skips it if `skip_initial_rules` is set to
    ///    `true`.
    /// 2. Greedily find an inclusion maximal clique `cliq`
    /// 3. Removes all nodes in `cliq` from the graph
    /// 4. Adds `cliq`.len() - 1 to the lower bound
    /// 5. Repeat from 1. until no more clique is found
    /// 6. The solution (that is found by the reduction rules) size is added to the lower bound.
    /// 7. Add the lower bound from `.graph.small_disjunct_cycles_heuristic()` to the lower bound. 
    ///
    /// Returns a lower bound.
    ///
    /// # Panics
    /// Panics if the first rule of `rule_priority` is not `Rule::SimpleRules`.
    pub fn lower_bound_clique_heuristic(&self, rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> usize {
        let mut clone_instance = self.clone();
        let mut lower_bound = 0;
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        assert!(!rule_priority.contains(&Rule::AdvancedPetal));
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            // Greedily find clique `cliq`
            let cliq = clone_instance.graph.greedy_max_clique();
            if cliq.len() > 1 {
                // TODO: remove cliques with higher degree first (especially when `cliq.len() ==
                // 2`.
                // Remove `cliq` from graph 
                clone_instance.graph.remove_nodes(cliq.clone());
                lower_bound += cliq.len() - 1;
            } else {
                // until no more clique is found.
                break
            }
        }
        // add solution size to both lower and upper.
        lower_bound += clone_instance.solution.len();
        // find bounds for left overs
        let instance = DFVSInstance::new(clone_instance.graph.clone(), None, None);
        // TODO: could be improved:
        let add_lower_sol = instance.disjunct_cycle_heuristic(rule_priority, true);
        lower_bound += add_lower_sol;
        return lower_bound
    }

    /// Lower bound heuristic that repeatedly applies simple reduction rules, finds a small cycle,
    /// removes that cycle, and adds 1 to the lower bound until no more cycle remains.
    /// 
    /// Returns a lower bound.
    pub fn disjunct_cycle_heuristic(&self, rule_priority: &Vec<Rule>, mut skip_initial_rules: bool) -> usize {
        let mut clone_instance = self.clone();
        let mut lower_bound = 0;
        // Simple rules are needed to check for circles. In addition it makes 
        // a lot of sense to have them in front of the `priority_list`.
        assert_eq!(rule_priority[0], Rule::SimpleRules);
        assert!(!rule_priority.contains(&Rule::AdvancedPetal));
        loop {
            // perform reductions:
            if !skip_initial_rules {
                clone_instance.exhaustive_reductions(rule_priority);
            } else {
                skip_initial_rules = false;
            }
            if clone_instance.graph.num_nodes() != 0 {
                let node = clone_instance.graph.nodes().next().expect("`clone_instance.graph` contains nodes");
                let mut cycle = clone_instance.graph.smallest_simple_cycle(node).expect("`clone_instance.graph` consists of cycles.");
                cycle.push(node);
                clone_instance.graph.remove_nodes(cycle.clone());
                lower_bound += 1;
            } else {
                break
            }
        }
        lower_bound += clone_instance.solution.len();
        lower_bound
    }

    // TODO: --------From here on we need to use finallize solution in different heursitics

    ///// Heuristic that does the following: 
    ///// 1. Apply some simple reduction rules
    ///// 2. Greedily find an inclusion maximal clique `cliq`
    ///// 3. Branch on clique: 
    /////     3.1. Removing and adding all nodes in `cliq` from the graph.
    /////     3.2. Removing and adding all but the node with the lowest weight of `weight_fn` given `weights`
    /////     `v` in `cliq` and contracting `v`.
    ///// 4. Repeat each branch from 1. until no more clique is found
    ///// 5. Add the solution from `left_over_heur` with `weight_fn` to the solution.
    ///// TODO: Make it possible to skip initial reduction. (like BST)
    /////
    ///// Returns the found solution, or `None` if execution was interrupted.
    /////
    ///// TODO: Takes ages, needs to operate on a few cliques only
    //pub fn advanced_clique_heuristic(
    //    &self, 
    //    left_over_heur: &dyn Fn(&DFVSInstance, &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, (f64, f64), &Vec<Rule>, bool) -> Option<FxHashSet<usize>>, 
    //    weight_fn: &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, 
    //    weights: (f64, f64), 
    //    rule_priority: &Vec<Rule>
    //    ) -> Result<Option<FxHashSet<usize>>, ProcessingError> {
    //    let mut clone_instance = self.clone();
    //    clone_instance._reset_changes();
    //    // Simple rules are needed to check for circles. In addition it makes 
    //    // a lot of sense to have them in front of the `priority_list`.
    //    assert_eq!(rule_priority[0], Rule::SimpleRules);
    //    clone_instance.advanced_clique_heuristic_recursion(left_over_heur, weight_fn, weights, rule_priority, 0)
    //}

    //fn advanced_clique_heuristic_recursion(
    //    &mut self, 
    //    left_over_heur: &dyn Fn(&DFVSInstance, &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, (f64, f64), &Vec<Rule>, bool) -> Option<FxHashSet<usize>>, 
    //    weight_fn: &dyn Fn(&Digraph, usize, (f64, f64)) -> Option<f64>, 
    //    weights: (f64, f64), 
    //    rule_priority: &Vec<Rule>,
    //    depth: usize
    //    ) -> Result<Option<FxHashSet<usize>>, ProcessingError> {
    //    self.start_new_changes();
    //    if self.check_interrupt() {
    //        return Ok(None)
    //    }
    //    // perform reductions:
    //    self.exhaustive_reductions(rule_priority);
    //    // Greedily find clique `cliq`
    //    let mut cliq = self.graph.greedy_max_clique();
    //    if cliq.len() > 4 && depth < 4 {
    //        self.start_new_changes();
    //        // Remove `cliq` from graph 
    //        self.add_all_to_solution(&cliq.clone()).expect("All those nodes exist.");
    //        let opt1;
    //        if let Some(mut opt1fu) = self.advanced_clique_heuristic_recursion(left_over_heur, weight_fn, weights, rule_priority, depth + 1)? {
    //            opt1fu.extend(cliq.clone());
    //            opt1 = opt1fu;
    //        } else {
    //            return Ok(None)
    //        }
    //        self.rebuild_section()?;
    //        let min_node = self.graph.get_min_weight_node_within(weight_fn, weights, &cliq).expect("|`cliq`| >= 1");
    //        cliq.remove(&min_node);
    //        self.add_all_to_solution(&cliq.clone()).expect("All those nodes exist.");
    //        self.graph.contract_node(min_node);
    //        let opt2;
    //        if let Some(mut opt2fu) = self.advanced_clique_heuristic_recursion(left_over_heur, weight_fn, weights, rule_priority, depth + 1)? {
    //            opt2fu.extend(cliq.clone());
    //            opt2 = opt2fu;
    //        } else {
    //            return Ok(None)
    //        }
    //        self.rebuild_section()?;
    //        if opt1.len() < opt2.len() {
    //            return Ok(Some(opt1))
    //        } else {
    //            return Ok(Some(opt2))
    //        }
    //    } else {
    //        // until no more clique is found.
    //        if self.graph.num_nodes() == 0 {
    //            self.rebuild_section()?;
    //            return Ok(Some(self.solution.clone()))
    //        }
    //        // find solution for left overs
    //        if let Some(mut upper_sol) = left_over_heur(&self, weight_fn, weights, rule_priority, true) {
    //            upper_sol.extend(self.solution.clone());
    //            self.rebuild_section()?;
    //            return Ok(Some(upper_sol))
    //        } 
    //        return Ok(None)
    //    }
    //}

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digraph::Digraph;
    use crate::dfvs_instance::DFVSInstance;
    use std::io::Cursor;

    #[test]
    fn weight_heuristic_is_deg_test() {
        let gr = Cursor::new("6 12 0\n2 6\n4 5\n2 5\n1 3\n4 6\n1 3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let instance = DFVSInstance::new(g, None, None);
        let sol = instance.top_down_weight_heuristic(&Digraph::lin_weight, (0f64,1f64), &vec![Rule::SimpleRules], false);
        assert_eq!(sol, vec![0usize, 2].into_iter().collect::<FxHashSet<usize>>());
    }

    #[test]
    fn small_cycle_heuristic_test() {
        let gr = Cursor::new("7 14 0\n2 3\n3 4\n4 5\n5 6\n6 7\n7 1\n1 2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let dvfsi = DFVSInstance::new(g, None, None);
        let lower_bound = dvfsi.disjunct_cycle_heuristic(&vec![Rule::SimpleRules], false);
        assert_eq!(lower_bound, 1);
    }

    #[test]
    fn local_search_test() {
        let gr = Cursor::new("5 9 0\n2 3\n1 3\n1 2 4\n5\n2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let dvfsi = DFVSInstance::new(g, None, None);
        let solution: FxHashSet<usize> = vec![0,1,3,4].into_iter().collect();
        assert!(dvfsi.local_search(&solution).is_some());
        let better_solution = dvfsi.exhaustive_local_search(&solution);
        assert_eq!(better_solution, 
                   Some(vec![0,1].into_iter().collect::<FxHashSet<usize>>()));
    }

    #[test]
    fn via_vertex_cover_test() {
        let gr = Cursor::new("7 28 0\n2 3 7 6\n1 3 4 7\n1 2 4 5\n2 3 5 6\n3 4 6 7\n\
                              4 5 7 1\n5 6 1 2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let dvfsi = DFVSInstance::new(g, None, None);
        if let VCQ::Exact(sol) = dvfsi.via_vertex_cover() {
            assert_eq!(sol.len(), 5);
        } else {
            assert!(false);
        }
        let gr = Cursor::new("5 12 0\n2 3 4 5\n1 4\n1 2\n1 5\n1 3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let dvfsi = DFVSInstance::new(g, None, None);
        if let VCQ::Bounds(bound,_) = dvfsi.via_vertex_cover() {
            assert_eq!(bound, 1);
        } else {
            assert!(false);
        }
    }

    // TODO: Are there missing tests? Please complete:
}
