//! 
//! This datastructure allows for many different alterations of the underlaying graph datastructure,
//! which are all saved in `reductions` and can be redone with `rebuild_section()` which redoes the 
//! last entries in `reductions` up to the index which is saved last in `register`.

use crate::digraph::Digraph;
use crate::reduction_rules::Rule;
use fxhash::FxHashSet;
use crate::cust_errors::{ImportError, ProcessingError};
use std::collections::HashSet;
use std::io::prelude::*;
use std::io;
use std::cmp::max;

/// All possible reductions that can be recorded.
#[derive(Debug, Clone, Eq, PartialEq)]
enum Reduction {
    /// Replaced node by edge
    Replace(usize, (usize, usize)), 
    /// Merged first node into second. Holds outlist of first node and introduced double edges.
    MergeBack(usize, usize, FxHashSet<usize>, FxHashSet<usize>), 
    /// Merged first node into second. Holds inlist of first node and introduced double edges.
    MergeForward(usize, usize, FxHashSet<usize>, FxHashSet<usize>),
    /// Merge the first node into the second. Holds the in- and outlist of the first, all the in- 
    /// and out neighbors both nodes have in common and the edges between the first and second node.
    CompleteMerge(usize, usize, (FxHashSet<usize>, FxHashSet<usize>), 
                  (FxHashSet<usize>, FxHashSet<usize>), HashSet<(usize, usize)>), 
    /// Contracts node. Also holds the incoming neigbors, the outgoing neigbors and the introduced 
    /// double edges.
    Contract(usize, FxHashSet<usize>, FxHashSet<usize>, HashSet<(usize, usize)>), 
    /// Added node to the solution and removed it. Additionally the incoming and outgoing neighbors
    /// are saved.
    AddedNode(usize, (FxHashSet<usize>, FxHashSet<usize>)),
    /// Removed node. Additionally the incoming and outgoing neighbors are saved.
    RemovedNode(usize,(FxHashSet<usize>, FxHashSet<usize>)), 
    /// Removed edge.
    RemovedEdge((usize, usize)), 
}

#[derive(Debug, Default, Clone, Eq, PartialEq)]
pub struct DFVSInstance {
    pub graph: Digraph,
    pub solution: FxHashSet<usize>,
    pub upper_bound: Option<usize>,
    pub lower_bound: Option<usize>,
    pub current_best: Option<FxHashSet<usize>>,
    pub merge_nodes: Vec<(usize, (usize, usize))>,
    /// Keeps the indcies of the first reductions in the different reduction steps.
    register: Vec<usize>,
    reductions: Vec<Reduction>,
}

impl DFVSInstance {

    /// Returns a new `DFVSInstance`
    ///
    /// # Arguments
    ///
    /// * `upper_bound` - If not specified (= None) the upper bound is computed by
    /// `DirectedFeedbackVertexSetInstance::high_degree_heuristic`.
    pub fn new(graph: Digraph, upper_bound: Option<usize>, lower_bound: Option<usize>) -> Self {
        DFVSInstance {
            graph,
            solution: FxHashSet::default(),
            upper_bound,
            lower_bound,
            current_best: None,
            merge_nodes: Vec::new(),
            register: vec![0],
            reductions: Vec::new(),
        }
    }

    /// Returns the effective upper bound of the current instance: 
    /// `self.upper_bound - self.solution.len()`
    pub fn effective_upper_bound(&self) -> Option<usize> {
        self.upper_bound.map(|upper| upper - self.solution.len())
    }

    /// Returns the effective lower bound of the current instance: 
    /// `self.lower_bound - self.solution.len()`
    pub fn effective_lower_bound(&self) -> Option<usize> {
        self.lower_bound.map(|lower| lower - self.solution.len())
    }
    
    /// Updates the old upper bound with `new_upper`. Or adds `new_upper` as upper bound if none
    /// exists.
    pub fn update_upper_bound(&mut self, new_upper: usize) {
        if let Some(upper) = self.upper_bound {
            if upper > new_upper {
                self.upper_bound = Some(new_upper);
            }
        } else {
            self.upper_bound = Some(new_upper);
        }
    }

    /// Updates the old lower bound with `new_lower`. Or adds `new_lowet` as lower bound if none
    /// exists.
    pub fn update_lower_bound(&mut self, new_lower: usize) -> bool {
        if let Some(lower) = self.lower_bound {
            if lower < new_lower {
                self.lower_bound = Some(new_lower);
                return true
            }
        } else {
            self.lower_bound = Some(new_lower);
            return true
        }
        false
    }

    /// Sets `set` to `self.current_best` and `self.upper_bound` to `set.len()` if `self.upper_bound` is 
    /// greater or equal to `set.len()` or `self.upper_bound` is not yet set. 
    /// Does nothing otherwise. 
    ///
    /// Returns `true` if `self.current_best` was updated and `false` otherwise.
    pub fn set_current_best(&mut self, set: &FxHashSet<usize>) -> bool {
        if let Some(upper) = self.upper_bound {
            if upper >= set.len() {
                self.current_best = Some(set.clone());
                self.upper_bound = Some(set.len());
                return true
            }
        } else {
            self.current_best = Some(set.clone());
            self.upper_bound = Some(set.len());
            return true
        }
        false
    }

    /// Computes a fast upper bound and a fitting solution. If the bound is better than the current
    /// best, sets a new upper bound.
    pub fn compute_and_set_fast_upper(&mut self, skip_initial_rules: bool) {
        let upper = self.get_fast_upper(skip_initial_rules);
        self.set_current_best(&upper);
    }

    /// Computes and sets the best currently available upper- and lower bounds and the best current
    /// solution.
    ///
    /// To access the best current solution use `self.current_best`.
    pub fn compute_and_set_upper_lower(&mut self, skip_initial_rules: bool) {
        let (lower, _, set) = self.get_best_bounds(skip_initial_rules);
        self.set_current_best(&set);
        self.update_lower_bound(lower);
    }

    /// Computes and sets a lower bound, if the computed lower bound is better then the old one. 
    pub fn compute_and_set_lower(&mut self, skip_initial_rules: bool) {
        let lower = self.get_some_lower(skip_initial_rules);
        self.update_lower_bound(lower); 
    }

    pub fn _reset_changes(&mut self) {
        self.register = vec![0];
        self.reductions = Vec::new();
    }

    /// This fuction is used to fix certain segments that can be recovered indipendently.
    pub fn start_new_changes(&mut self) {
        self.register.push(self.reductions.len());
    }

    /// Rebuilds the last reduction segment.
    /// Returns an error if the recovery failed.
    pub fn rebuild_section(&mut self) -> Result<(), ProcessingError> {
        if self.register.is_empty() {
            return Err(ProcessingError::RebuildError)
        }
        let up_to = self.register.pop().expect("We checked register.");
        while self.reductions.len() > up_to {
            match self.reductions.pop().expect("This can not be empty.") {
                Reduction::Replace(node, edge) => {
                    if !self.graph.remove_edge(&edge) || 
                        !self.graph.reinsert_node(
                            node, 
                            vec![edge.0].into_iter().collect(), 
                            vec![edge.1].into_iter().collect()
                            ) {
                            return Err(ProcessingError::RebuildError);
                    }
                },
                Reduction::RemovedEdge(edge) => self.graph.add_edge(edge),
                Reduction::RemovedNode(node, (ins, outs)) => {
                    if !self.graph.reinsert_node(node, ins, outs) {
                        return Err(ProcessingError::RebuildError)
                    }
                },
                Reduction::AddedNode(node, (ins, outs)) => {
                    if !self.graph.reinsert_node(node, ins, outs) {
                        return Err(ProcessingError::RebuildError)
                    }
                    self.solution.remove(&node);
                },
                Reduction::MergeBack(node, into, outs, doubles) => {
                    // reinsert node
                    if !self.graph.reinsert_node(node, vec![into].into_iter().collect(), outs.clone()){
                        return Err(ProcessingError::RebuildError)
                    }
                    // collect where node needs to be removed 
                    let remove_from_ins = outs.difference(&doubles);
                    // remove node from lists 
                    self.graph.remove_edges(remove_from_ins.map(|trg| (into, *trg)));
                },
                Reduction::MergeForward(node, into, ins, doubles) => {
                    // reinsert node
                    if !self.graph.reinsert_node(node, ins.clone(), vec![into].into_iter().collect()){
                        return Err(ProcessingError::RebuildError)
                    }
                    // collect where node needs to be removed 
                    let remove_from_outs = ins.difference(&doubles);
                    // remove node from lists 
                    self.graph.remove_edges(remove_from_outs.map(|src| (*src, into)));
                },
                Reduction::CompleteMerge(node, into, singles, doubles, betweens) => {
                    // reinsert node 
                    if !self.graph.reinsert_node(node, singles.0.clone(), singles.1.clone()) {
                        return Err(ProcessingError::RebuildError)
                    }
                    for edge in betweens {
                        self.graph.add_edge(edge);
                    }
                    // collect nodes that need to be removed from `into`
                    let remove_from_outs = singles.0.difference(&doubles.0);
                    let remove_from_ins = singles.1.difference(&doubles.1);
                    self.graph.remove_edges(remove_from_outs.map(|src| (*src, into)));
                    self.graph.remove_edges(remove_from_ins.map(|trg| (into, *trg)));
                    let placeholder = self.merge_nodes.len() + self.graph.num_reserved_nodes() - 1;
                    self.merge_nodes.pop();
                    self.solution.remove(&placeholder);
                }
                Reduction::Contract(node, ins, outs, doubles) => {
                    if !self.graph.reinsert_node(node,ins.clone(), outs.clone()) {
                        return Err(ProcessingError::RebuildError)
                    }
                    // remove all but doubles
                    self.graph.remove_edges(ins.iter()
                                            .flat_map(|src| {
                                                outs.iter().map(move |trg| (*src, *trg))
                                            })
                                            .collect::<HashSet<_>>()
                                            .difference(&doubles)
                                            .copied());
                }
            }
        }
        if self.register.is_empty() {
            self.register.push(0);
        }
        Ok(())
    }

    /// TODO Needs to be incorperated into the new redo changes
    /// Reverses the changes completely.
    pub fn rebuild_complete(&mut self) -> Result<(), ProcessingError> {
        if self.register.is_empty() {
            return Err(ProcessingError::RebuildError)
        }
        for _ in 0..self.register.len() {
            self.rebuild_section()?;
        }
        Ok(())
    }

    /// Returns current state of the registers.
    pub fn check_registers(&self) -> (usize, Vec<usize>) {
        (self.reductions.len(), self.register.clone())
    }

    /// Removes `node` from the graph.
    /// Throws an error is `node` was not in the graph.
    pub fn remove_node(&mut self, node: usize) -> Result<(), ProcessingError>{
        if let Some((ins, outs)) = self.graph.remove_node(node) {
            self.reductions.push(Reduction::RemovedNode(node, (ins, outs)));
            return Ok(())
        }
        return Err(ProcessingError::InvalidParameter("Given node was not contained in the graph.".to_owned()))
    }

    /// Removes `node` from the graph.
    /// Returns effected nodes.
    /// Throws an error if `node` was not in the graph.
    pub fn remove_node_return_effected(&mut self, node: usize) -> 
        Result<FxHashSet<usize>, ProcessingError>{
        let mut effected = FxHashSet::default();
        if let Some((ins, outs)) = self.graph.remove_node(node) {
            effected.extend(ins.iter());
            effected.extend(outs.iter());
            self.reductions.push(Reduction::RemovedNode(node, (ins, outs)));
            return Ok(effected)
        }
        return Err(ProcessingError::InvalidParameter("Given node was not contained in the graph.".to_owned()))
    }

    /// Removes all nodes in `node_set` from `self.graph`. 
    /// Returns `Ok` and records the alteration if all nodes were added, returns a `ProcessingError` otherwise.
    pub fn remove_nodes(&mut self, node_set: &FxHashSet<usize>) -> Result<(), ProcessingError> {
        for node in node_set {
            self.remove_node(*node)?;
        }
        Ok(())
    }

    /// Removes all nodes in `node_set` from `self.graph`. 
    /// Returns the effected nodes, that remain in `self.graph` and records the alteration if all nodes were added, returns a `ProcessingError` otherwise.
    pub fn remove_nodes_return_effected(&mut self, node_set: &FxHashSet<usize>) -> Result<FxHashSet<usize>, ProcessingError> {
        let mut effected = FxHashSet::default();
        for node in node_set {
            if let Some((ins, outs)) = self.graph.remove_node(*node) {
                effected.extend(ins.iter());
                effected.extend(outs.iter());
                self.reductions.push(Reduction::RemovedNode(*node, (ins, outs)));
            } else {
                return Err(ProcessingError::InvalidParameter("Given node set was not completely contained in the graph.".to_owned()))
            }
        }
        Ok(effected.into_iter().filter(|n| !node_set.contains(n)).collect())
    }

    /// Adds `node` to the solution and removes it from the graph.
    /// Throws an error is `node` was not in the graph.
    pub fn add_to_solution(&mut self, node: usize) -> Result<(), ProcessingError>{
        if let Some((ins, outs)) = self.graph.remove_node(node) {
            self.solution.insert(node);
            self.reductions.push(Reduction::AddedNode(node, (ins, outs)));
            return Ok(())
        }
        return Err(ProcessingError::InvalidParameter("Given node was not contained in the graph.".to_owned()))
    }

    /// Adds `node` to the solution and removes it from the graph. 
    /// Returns effected nodes.
    /// Throws an error is `node` was not in the graph.
    pub fn add_to_solution_return_effected(&mut self, node: usize) -> 
        Result<FxHashSet<usize>, ProcessingError>{
        let mut effected = FxHashSet::default(); 
        if let Some((ins, outs)) = self.graph.remove_node(node) {
            self.solution.insert(node);
            effected.extend(ins.iter());
            effected.extend(outs.iter());
            self.reductions.push(Reduction::AddedNode(node, (ins, outs)));
            return Ok(effected)
        }
        return Err(ProcessingError::InvalidParameter("Given node was not contained in the graph.".to_owned()))
    }

    /// Adds all nodes in `node_set` to `self.solution` and removes them from `self.graph`. 
    /// Returns `Ok` and records the alteration if all nodes were added, returns a `ProcessingError` otherwise.
    pub fn add_all_to_solution(&mut self, node_set: &FxHashSet<usize>) -> Result<(), ProcessingError> {
        for node in node_set {
            if let Some((ins, outs)) = self.graph.remove_node(*node) {
                self.reductions.push(Reduction::AddedNode(*node, (ins, outs)));
                self.solution.insert(*node);
            } else {
                return Err(ProcessingError::InvalidParameter("Given node set was not completely contained in the graph.".to_owned()))
            }
        }
        Ok(())
    }

    /// Adds all nodes in `node_set` to `self.solution` and removes them from `self.graph`. 
    /// Returns all effected nodes still in `self.graph` and records the alteration if all nodes were added, returns a `ProcessingError` otherwise.
    pub fn add_all_to_solution_return_effected(&mut self, node_set: &FxHashSet<usize>) -> Result<FxHashSet<usize>, ProcessingError> {
        let mut effected = FxHashSet::default();
        for node in node_set {
            if let Some((ins, outs)) = self.graph.remove_node(*node) {
                effected.extend(ins.iter());
                effected.extend(outs.iter());
                self.reductions.push(Reduction::AddedNode(*node, (ins, outs)));
                self.solution.insert(*node);
            } else {
                return Err(ProcessingError::InvalidParameter("Given node set was not completely contained in the graph.".to_owned()))
            }
        }
        Ok(effected.into_iter().filter(|n| !node_set.contains(n)).collect())
    }

    /// Removes `edge` and adds `node` with `edge.0` as an incoming neighbor and `edge.1` as an
    /// outgoing neighbor and records the alteration.
    /// Throws an error is `node` was not in the graph.
    pub fn replace(&mut self, node: usize, edge: (usize,  usize)) -> Result<(), ProcessingError> {
        if self.graph.remove_node(node).is_some() {
            // check if edge is replaced (double edge)
            if self.graph.edge_exists(&edge) {
                self.reductions
                    .push(
                        Reduction::RemovedNode(node, 
                                               (vec![edge.0]
                                                .into_iter()
                                                .collect::<FxHashSet<usize>>(),
                                               vec![edge.1]
                                               .into_iter()
                                               .collect::<FxHashSet<usize>>()
                                               )));
            } else {
                self.graph.add_edge(edge);
                self.reductions.push(Reduction::Replace(node, edge));
            }
            return Ok(())
        }
        return Err(ProcessingError::InvalidParameter("Given node was not contained in the graph.".to_owned()))
    }

    /// Merges `node` with it's only incoming neighbor `into`.
    /// Throws an error is `node` was not in the graph.
    pub fn merge_into_back(&mut self, node: usize, into: usize) -> Result<(), ProcessingError> {
        if let Some((outs, doubles)) = self.graph.merge_into_back(node, into) {
            self.reductions.push(Reduction::MergeBack(node, into, outs, doubles));
            return Ok(())
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }

    /// Merges `node` with it's only outgoing neighbor `into`.
    /// Throws an error is `node` was not in the graph.
    pub fn merge_into_front(&mut self, node: usize, into: usize) -> Result<(), ProcessingError> {
        if let Some((ins, doubles)) = self.graph.merge_into_front(node, into) {
            self.reductions.push(Reduction::MergeForward(node, into, ins, doubles));
            return Ok(())
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }

    /// Merges `node` with it's only incoming neighbor `into`.
    /// Returns the outgoing neighbors.
    /// Throws an error is `node` was not in the graph.
    pub fn merge_into_back_return_front(&mut self, node: usize, into: usize) -> Result<FxHashSet<usize>, ProcessingError> {
        if let Some((outs, doubles)) = self.graph.merge_into_back(node, into) {
            self.reductions.push(Reduction::MergeBack(node, into, outs.clone(), doubles));
            return Ok(outs)
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }

    /// Merges `node` with it's only outgoing neighbor `into`.
    /// Returns the incoming neighbors.
    /// Throws an error is `node` was not in the graph.
    pub fn merge_into_front_return_back(&mut self, node: usize, into: usize) -> Result<FxHashSet<usize>, ProcessingError> {
        if let Some((ins, doubles)) = self.graph.merge_into_front(node, into) {
            self.reductions.push(Reduction::MergeForward(node, into, ins.clone(), doubles));
            return Ok(ins)
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }


    /// Contracts `link` and merges `neighbors[0]` into `neighbors[1]`. Push
    /// `self.graph.num_reserved_nodes() + self.merge_nodes.len()` into the solution as a placeholder, and
    /// add information to `self.merge_nodes` to figure how the placeholder will be converted.
    pub fn contract_link_node(&mut self, link: usize, neighbors: &[usize; 2]) -> Result<(), ProcessingError> {
        self.remove_node(link)?;
        if let Some((single_neighs, doubles, betweens)) = self.graph.complete_merge(neighbors[0], neighbors[1]) {
            // remove a node from the solution and pop placeholder.
            self.reductions.push(Reduction::CompleteMerge(neighbors[0], neighbors[1], single_neighs, doubles, betweens));
            let id = self.merge_nodes.len() + self.graph.num_reserved_nodes();
            self.solution.insert(id);
            self.merge_nodes.push((neighbors[1], (neighbors[0], link)));
            return Ok(())
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }

    /// Contracts `link` and merges `neighbors[0]` into `neighbors[1]`. Push
    /// `self.graph.num_reserved_nodes() + self.merge_nodes.len()` into the solution as a placeholder, and
    /// add information to `self.merge_nodes` to figure how the placeholder will be converted.
    /// Returns the effected node that still remain in `self.graph`.
    pub fn contract_link_node_return_effected(&mut self, link: usize, neighbors: &[usize; 2]) -> Result<(), ProcessingError> {
        self.remove_node(link)?;
        if let Some((single_neighs, doubles, betweens)) = self.graph.complete_merge(neighbors[0], neighbors[1]) {
            // remove a node from the solution and pop placeholder.
            self.reductions.push(Reduction::CompleteMerge(neighbors[0], neighbors[1], single_neighs, doubles, betweens));
            let id = self.merge_nodes.len() + self.graph.num_reserved_nodes();
            self.solution.insert(id);
            self.merge_nodes.push((neighbors[1], (neighbors[0], link)));
            return Ok(())
        } 
        Err(ProcessingError::InvalidParameter("Given nodes were not suited for this merge operation.".to_owned()))
    }

    /// Removes `edge` from the graph.
    /// Records the alteration.
    pub fn remove_edge(&mut self, edge: &(usize, usize)) -> Result<(), ProcessingError> {
        if self.graph.remove_edge(edge) {
            self.reductions.push(Reduction::RemovedEdge(*edge));
            return Ok(())
        }
        Err(ProcessingError::InvalidParameter("Given edge did not exist.".to_owned()))
    }
    
    /// Removes directed edges of the form `(src, trg)` from the graph
    /// Records the alterations
    ///
    /// # Panics
    /// Panics if a node index is invalid.
    pub fn remove_edges<I: IntoIterator<Item=(usize, usize)>>(&mut self, edges: I) -> Result<(), ProcessingError> {
        for edge in edges {
            self.remove_edge(&edge)?;
        }
        Ok(())
    }

    /// Contracts `node` adding edges between the incoming and the outgoing neighbors of `node`. 
    /// Returns `true` if `node` was contracted, and `false` if not.
    /// Records the alteration.
    ///
    /// # Panics 
    /// Panics if `node` is out of bounds.
    pub fn contract_node(&mut self, node: usize) -> Result<(), ProcessingError> {
        if let Some((ins, outs, doubles)) = self.graph.contract_node(node) {
            self.reductions.push(Reduction::Contract(node, ins, outs, doubles));
            Ok(())
        } else {
            Err(ProcessingError::InvalidParameter("Given node did not exist.".to_owned()))
        }
    }

    /// Replaces placeholder from `self.contract_link_node()` with the actual intended nodes.
    /// Can't currently be reversed.
    /// TODO: Make it not fuck up the reverse.
    pub fn finallize_solution(&mut self) {
        while !self.merge_nodes.is_empty() {
            let id = self.graph.num_reserved_nodes() + self.merge_nodes.len() - 1;
            self.solution.remove(&id);
            let rule = self.merge_nodes.pop().expect("`merge_copy` is not empty");
            if self.solution.contains(&rule.0) {
                self.solution.insert(rule.1.0);
            } else {
                self.solution.insert(rule.1.1);
            }
        }
    }

    /// Finalizes the solution as in `self.finallize_solution()` without actually changing the
    /// solution, so the actual instance can still be redone.
    pub fn finallize_solution_temp(&mut self) -> FxHashSet<usize> {
        let mut merge_copy = self.merge_nodes.clone();
        let mut solution_copy = self.solution.clone();
        while !merge_copy.is_empty() {
            let id = self.graph.num_reserved_nodes() + merge_copy.len() - 1;
            solution_copy.remove(&id);
            let rule = merge_copy.pop().expect("`merge_copy` is not empty");
            if solution_copy.contains(&rule.0) {
                solution_copy.insert(rule.1.0);
            } else {
                solution_copy.insert(rule.1.1);
            }
        }
        solution_copy
    }

    /// Finalizes a given solution `sol` as in `self.finallize_solution()` without actually changing
    /// `self.solution`, so the actual instance can still be redone.
    pub fn finallize_given_solution_temp(&self, sol: &FxHashSet<usize>) -> FxHashSet<usize> {
        let mut merge_copy = self.merge_nodes.clone();
        let mut solution_copy = sol.clone();
        while !merge_copy.is_empty() {
            let id = self.graph.num_reserved_nodes() + merge_copy.len() - 1;
            solution_copy.remove(&id);
            let rule = merge_copy.pop().expect("`merge_copy` is not empty");
            if solution_copy.contains(&rule.0) {
                solution_copy.insert(rule.1.0);
            } else {
                solution_copy.insert(rule.1.1);
            }
        }
        solution_copy
    }

    /// Finds the best upper- and lower bound under the currently implemented heuristics (or
    /// approximations).
    ///
    /// Returns the lower bound, the upper bound and the solution of the best upper bound.
    ///
    /// TODO: make more customizable if we have more time. 
    /// TODO: Update rules
    pub fn get_best_bounds(&self, skip_initial_rules: bool) -> (usize, usize, FxHashSet<usize>) {
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        let rule_priority = &vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC];
        let bounds = self.clique_heuristic(&rule_priority, skip_initial_rules);
        upper.push(bounds.clone().2);
        lower.push(bounds.clone().0);
        upper.push(self.top_down_weight_heuristic(&Digraph::cai_weight, (0.2,0f64), &rule_priority, skip_initial_rules));
        // TODO: this might be a waste of time. Find better lower bounds:
        //lower.push(self.disjunct_cycle_heuristic(&rule_priority, skip_initial_rules));
        let mut best_upper = upper.iter().max_by_key(|ub| ub.len()).expect("There is at least one upper bound").clone();

        if let Some(better) = self.exhaustive_local_search(&best_upper) {
            best_upper = better;
        }
        let best_lower = lower.iter().min().expect("There is at least one lower bound.");
        return (*best_lower, best_upper.len(), best_upper) 
    }

    /// Returns a good upper bound for the given instance.
    pub fn get_good_upper(&self, skip_initial_rules: bool) -> FxHashSet<usize> {
        let rule_priority = &vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::Crown];
        let mut upper = self.top_down_weight_heuristic(&Digraph::cai_weight, (0.2,0f64), &rule_priority, skip_initial_rules);
        if let Some(better) = self.exhaustive_local_search(&upper) {
            upper = better;
        }
        return upper
    }

    /// Returns a upper bound for the given instance.
    pub fn get_fast_upper(&self, skip_initial_rules: bool) -> FxHashSet<usize> {
        let rule_priority = &vec![Rule::SimpleRules];
        let mut upper = self.top_down_weight_heuristic_only_local_simple(&Digraph::cai_weight, (0.2,0f64), skip_initial_rules);
        return upper
    }

    /// Finds the best lower bound under the currently implemented heuristics (or
    /// approximations).
    /// Applied heuristics only use simple rules.
    ///
    /// Returns the best lower bound found.
    pub fn get_some_lower(&self, skip_initial_rules: bool) -> usize {
        let (lower1, _, _) = self.clique_heuristic(&vec![Rule::SimpleRules], skip_initial_rules);
        let lower2 = self.disjunct_cycle_heuristic(&vec![Rule::SimpleRules], skip_initial_rules);
        return max(lower2, lower1)
    }

}

impl DFVSInstance {

    /// Reads the solution of a DFVS instance from a `BufRead` type.
    /// Returns a HashSet of nodes in the solution.
    pub fn read_solution<R: BufRead>(sol: R) -> Result<FxHashSet<usize>, ImportError> {
        sol.lines()
            .map(|line| {
                line.unwrap().parse::<usize>().or(Err(ImportError::InputMalformedError))
            }).collect()
    }

    /// Writes a solution to a `Write` type.
    pub fn write_solution<W: Write>(solution: &FxHashSet<usize>, mut out: W) -> Result<(), io::Error> { 
        for elem in solution {
            writeln!(out, "{}",elem + 1)?;
        }
        Ok(())
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digraph::Digraph;
    use std::io::Cursor;
    use fxhash::FxHashSet;
        
    // TODO: Test whether reduced and unreduced lead to the same result
    #[test]
    fn read_write_sol_test() {
        let gr = Cursor::new("5 13 0\n3 5\n3 4\n1 2 4 5\n1 3 5\n2 3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let instance = DFVSInstance::new(g, None, None);
        let sol = "2\n4";
        let sol_curs = Cursor::new(sol);
        let read_sol = DFVSInstance::read_solution(sol_curs);
        assert!(read_sol.is_ok());
        assert_eq!(read_sol.unwrap(), vec![2,4].into_iter().collect::<FxHashSet<usize>>());
        let stdout = io::stdout();
        let stdout = stdout.lock();
        assert!(DFVSInstance::write_solution(&instance.solution, stdout).is_ok());
    }

}
