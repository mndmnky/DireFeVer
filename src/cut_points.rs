//! 
//! This module implements different functions meant to find a set of nodes that cut the graph into
//! multiple strongly connected components.

use crate::digraph::Digraph;
use fxhash::{FxHashSet, FxHashMap};
use std::collections::HashMap;
use std::ops::AddAssign;
use std::vec;

impl Digraph {
    
    /// Finds a set of n nodes, such that their removal would split the graph into at least n+1
    /// strongly connected components.
    ///
    /// # Arguments 
    ///
    /// * `cut_size` - The amount of nodes tried in each step.
    pub fn find_split_set(&self, cut_size: usize) -> FxHashSet<usize> {
        // get all sccs 
        let subs = self.split_into_connected_components_alt();
        let mut cut_vertices: FxHashSet<usize> = FxHashSet::default();
        for mut sub in subs {
            if cut_size == 1 {
                if let Some(cut_vs) = sub.split_scc_iter() {
                    cut_vertices.extend(cut_vs.into_iter());
                }
            } else {
                if let Some(cut_vs) = sub.split_scc_recursively(cut_size) {
                    cut_vertices.extend(cut_vs.into_iter());
                }
            }
        }
        return cut_vertices;
    }

    /// Finds a set of n nodes, such that their removal would split the graph into strongly connected
    /// components.
    pub fn find_split_set_adv(&self) -> FxHashSet<usize> {
        // get all sccs 
        let subs = self.split_into_connected_components_alt();
        let mut cut_vertices: FxHashSet<usize> = FxHashSet::default();
        for mut sub in subs {
            if let Some(cut_vs) = sub.split_scc_adv() {
                cut_vertices.extend(cut_vs.into_iter());
            }
        }
        return cut_vertices;
    }

    /// TODO: Repeat and compare different levels.
    /// Finds a set of n nodes, such that their removal would split the graph into strongly connected
    /// components.
    ///
    /// # Arguments 
    ///
    /// * `cut_size` - The amount of nodes tried in each step.
    pub fn find_split_set_by_edge_contraction(&self, cut_size: usize) -> Option<(FxHashSet<usize>, usize)> {
        let mut gc = self.clone();
        let mut split_set = FxHashSet::default();
        let mut scc_count = 0;
        loop {
            if gc.num_edges() == 0 {
                break
            }
            let size = gc.num_nodes();
            let mut reps: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
            // Contract n/2 edges
            let mut gc_nodes: FxHashSet<usize> = gc.nodes().collect();
            while gc_nodes.len() > (size as f64 / 2.0).floor() as usize {
                // Find edge to contract 
                let mut low = gc.get_min_degree_node_within(&gc_nodes).expect("`gc` still has nodes");
                if let Some(mut lowest) = gc.get_min_degree_node_within(&gc.neighbors_in(low, &gc_nodes).expect("`low` exists.")) {
                    if !gc.has_edge((low, lowest)) {
                        let swap = low;
                        low = lowest;
                        lowest = swap;
                    }
                    let ev_lowest_node = if let Some(lowest_nodes) = reps.get(&lowest) {
                        Some(lowest_nodes.clone())
                    } else {
                        None
                    };
                    let low_rep = reps.entry(low).or_insert(Vec::new());
                    low_rep.push(lowest);
                    if let Some(lowest_reps) = ev_lowest_node {
                        low_rep.extend(lowest_reps.into_iter());
                    }
                    gc.contract_edge((low, lowest));
                    gc_nodes.remove(&low);
                    gc_nodes.remove(&lowest);
                } else {
                    break
                }
            }
            // Find bridges between sccs
            if let Some(cuts) = gc.split_scc_recursively_on_contracted(&reps, &mut scc_count, cut_size) {
                for cut in cuts {
                    split_set.insert(cut);
                    if reps.contains_key(&cut) {
                        split_set.extend(reps.get(&cut).expect("contains key"));
                    }
                }
                break
            }
        }
        if scc_count > 1 {
            return Some((split_set, scc_count))
        }
        None 
    }

    /// Looks for a node that splits `self` into two or more strongly connected components and
    /// repeats iteratively until no more node can be found. Returns the set of cut nodes, or
    /// `None` if none could have been found.
    fn split_scc_iter(&mut self) -> Option<FxHashSet<usize>> {
        let mut cut_vertices = FxHashSet::default();
        let mut parts = vec![self.clone()];
        while !parts.is_empty() {
            let part = parts.pop().expect("is not empty");
            let mut nodes: Vec<usize> = self.nodes().collect();
            nodes.sort_unstable_by_key(|n| self.degree(*n));
            while !nodes.is_empty() {
                let tn = nodes.pop().expect("`nodes` is not empty");
                let mut clone = part.clone();
                clone.remove_node(tn);
                if let Some(sccs) = clone.reduce_to_sccs() {
                    if sccs.len() > 1 {
                        cut_vertices.insert(tn);
                        parts.extend(clone.split_into_connected_components_alt().into_iter());
                        break;
                    }
                }
            }
        }
        if cut_vertices.len() == 0 {
            return None
        }
        return Some(cut_vertices)
    }

    /// Looks for a node that splits `self` into two or more strongly connected components and
    /// repeats recursively until no more node can be found. Returns the set of cut nodes, or
    /// `None` if none could have been found.
    ///
    /// # Arguments 
    ///
    /// * `cut_size` - The amount of nodes tried in each step.
    fn split_scc_recursively(&mut self, cut_size: usize) -> Option<FxHashSet<usize>> {
        let mut nodes: Vec<usize> = self.nodes().collect();
        nodes.sort_unstable_by_key(|n| self.degree(*n));
        let mut clone = self.clone();
        let mut cut_vertices = FxHashSet::default();
        for _ in 0..cut_size {
            if nodes.len() > 0 {
                let tn = nodes.pop().expect("`nodes` is not empty");
                clone.remove_node(tn);
                cut_vertices.insert(tn);
            }
        }
        if let Some(sccs) = clone.reduce_to_sccs() {
            if sccs.len() > 1 {
                let subs = clone.split_into_connected_components_alt();
                for mut sub in subs {
                    if let Some(cut_vs) = sub.split_scc_recursively(cut_size){
                        cut_vertices.extend(cut_vs.into_iter());
                    }
                }
                return Some(cut_vertices);
            }
        }
        None
    }

    /// Looks for a node that splits `self` into two or more strongly connected components and
    /// repeats recursively until no more node can be found. Returns the set of cut nodes, or
    /// `None` if none could have been found.
    ///
    /// # Arguments 
    ///
    /// * `cut_size` - The amount of nodes tried in each step.
    fn split_scc_recursively_on_contracted(&mut self, contraction_map: &FxHashMap<usize, Vec<usize>>, scc_amn: &mut usize, cut_size: usize) -> Option<FxHashSet<usize>> {
        if cut_size == 1 {
            let mut nodes: Vec<usize> = self.nodes().collect();
            nodes.sort_unstable_by_key(|n| self.degree(*n));
            nodes.reverse();
            nodes.sort_by_key(|n| {
                if contraction_map.contains_key(n) {
                    contraction_map.get(n).expect("contained").len()
                } else {
                    1
                }
            });
            nodes.reverse();
            while !nodes.is_empty() {
                let tn = nodes.pop().expect("`nodes` is not empty");
                let mut clone = self.clone();
                clone.remove_node(tn);
                if let Some(sccs) = clone.reduce_to_sccs_allow_contracted(contraction_map) {
                    let mut cut_vertices: FxHashSet<usize> = vec![tn].into_iter().collect();
                    if sccs.len() > 1 {
                        let subs = clone.split_into_connected_components_alt();
                        for mut sub in subs {
                            if let Some(cut_vs) = sub.split_scc_recursively_on_contracted(contraction_map, scc_amn, cut_size){
                                cut_vertices.extend(cut_vs.into_iter());
                            }
                        }
                        scc_amn.add_assign(sccs.len());
                        return Some(cut_vertices);
                    }
                }
            }
            return None
        } else {
            let mut nodes: Vec<usize> = self.nodes().collect();
            nodes.sort_unstable_by_key(|n| self.degree(*n));
            nodes.reverse();
            nodes.sort_by_key(|n| {
                if contraction_map.contains_key(n) {
                    contraction_map.get(n).expect("contained").len()
                } else {
                    1
                }
            });
            nodes.reverse();
            let mut clone = self.clone();
            let mut cut_vertices = FxHashSet::default();
            for _ in 0..cut_size {
                if nodes.len() > 0 {
                    let tn = nodes.pop().expect("`nodes` is not empty");
                    clone.remove_node(tn);
                    cut_vertices.insert(tn);
                }
            }
            if let Some(sccs) = clone.reduce_to_sccs_allow_contracted(contraction_map) {
                if sccs.len() > 1 {
                    let subs = clone.split_into_connected_components_alt();
                    for mut sub in subs {
                        if let Some(cut_vs) = sub.split_scc_recursively_on_contracted(contraction_map, scc_amn, cut_size){
                            cut_vertices.extend(cut_vs.into_iter());
                        }
                    }
                    scc_amn.add_assign(sccs.len());
                    return Some(cut_vertices);
                }
            }
            return None
        }
    }

    /// Looks for a node that splits `self` into two or more strongly connected components and
    /// repeats recursively until no more node can be found. Returns the set of cut nodes, or
    /// `None` if none could have been found.
    ///
    /// Allows for single nodes sccs if they have a loop.
    fn split_scc_recursively_allow_loops(&mut self, min_split: usize, running_count: &mut usize) -> Option<FxHashSet<usize>> {
        let mut nodes: Vec<usize> = self.nodes().collect();
        nodes.sort_unstable_by_key(|n| self.degree(*n));
        while !nodes.is_empty() {
            let tn = nodes.pop().expect("`nodes` is not empty");
            let mut clone = self.clone();
            clone.remove_node(tn);
            if let Some(sccs) = clone.reduce_to_sccs_allow_loops() {
                let mut cut_vertices: FxHashSet<usize> = vec![tn].into_iter().collect();
                if sccs.len() > min_split {
                    let subs = clone.split_into_connected_components_alt();
                    running_count.add_assign(sccs.len());
                    for mut sub in subs {
                        if let Some(cut_vs) = sub.split_scc_recursively_allow_loops(min_split, running_count){
                            cut_vertices.extend(cut_vs.into_iter());
                        }
                    }
                    return Some(cut_vertices);
                }
            }
        }
        None
    }

    /// Tries to split a graph into sccs with the help of `self.split_scc_recursively()`. If non
    /// cut vertices could been found, merge as much node disjoint 4 cycles as we can found into
    /// nodes, and try to split then.
    fn split_scc_adv(&mut self) -> Option<FxHashSet<usize>> {
        // find cuts as in the old version
        if let Some(cut_vs) = self.split_scc_recursively(1) {
            return Some(cut_vs);
        }
        // if none where found:
        // find indipendent cycles of size 4 
        let cycles = self.find_disjoint_cycles_of_at_most_size(4);
        let pairs = cycles.into_iter().map(|vec4| (vec4[0], vec4)).collect::<HashMap<usize, Vec<usize>>>();
        for (rep, cycle) in &pairs {
            // merge with loops
            self.big_merge_allow_loops(*rep, cycle[1..cycle.len()].into()).expect("all nodes in the cycles exist");
        }
        // find cuts 
        let mut scc_count = 0;
        let mut final_cuts = FxHashSet::default();
        if let Some(cut_vs) = self.split_scc_recursively_allow_loops(4, &mut scc_count) {
            for cut_v in cut_vs {
                if let Some(cycle) = pairs.get(&cut_v) {
                    final_cuts.extend(cycle.into_iter());
                }
            }
            eprintln!("ratio (cut_vs:sccs) {}:{}",final_cuts.len(), scc_count);
            return Some(final_cuts)
        }
        None
    }

}
