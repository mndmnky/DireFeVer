//! Directed, dynamic graph datastructure.
//! The main fields are the& two adjacency lists `in_list` and `out_list` which respectively hold
//! all incoming and outgoing neighbors of a node in a `FxHashSet`. 

use std::io::prelude::*;
use std::io;
use std::collections::{HashMap, HashSet, VecDeque};
use std::cmp::min;
use crate::cust_errors::{ImportError, ProcessingError};
use crate::other_ds::NodeSet;
use fxhash::{FxHashMap, FxHashSet};


/// Helper marker for the iterative approach to find all strongly connected components.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Im {
    Itm(usize),
    Marker(usize),
}

/// The graph datastructure
#[derive(Debug, Default, Clone, Eq, PartialEq)]
pub struct Digraph {
    in_list: Vec<Option<FxHashSet<usize>>>,
    out_list: Vec<Option<FxHashSet<usize>>>,
}

impl Digraph {

    /// Returns a hashset of the neighborhood of outgoing nodes of `node` if `node` was not already deleted.
    pub fn out_neighbors(&self, node: usize) -> &Option<FxHashSet<usize>> {
        &self.out_list[node]
    }

    /// Returns a hashset of the neighborhood of incoming nodes of `node` if `node` was not already deleted.
    pub fn in_neighbors(&self, node: usize) -> &Option<FxHashSet<usize>> {
        &self.in_list[node]
    }

    /// Returns a hashset of the neighborhood of `node` if `node` was not already deleted.
    pub fn neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let out_set = &self.out_list[node].as_ref().expect("`in_list` exists");
            Some(out_set.union(in_set).copied().collect::<FxHashSet<usize>>())
        } else {
            None
        }
    }

    /// Returns a hashset of the strongly connected neighborhood of `node`, or None if `node` was
    /// already removed.
    pub fn strong_neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let out_set = &self.out_list[node].as_ref().expect("`node` exists");
            let strong_neighbors = in_set.intersection(out_set).copied().collect::<FxHashSet<_>>();
            Some(strong_neighbors)
        } else {
            None
        }
    }

    /// Returns a hashset of the strongly connected closed neighborhood of `node`, or None if `node` was
    /// already removed.
    pub fn closed_strong_neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let out_set = &self.out_list[node].as_ref().expect("`node` exists");
            let mut strong_neighbors = in_set.intersection(out_set).copied().collect::<FxHashSet<_>>();
            strong_neighbors.insert(node);
            Some(strong_neighbors)
        } else {
            None
        }
    }

    /// Returns the set of weak incoming neighbors of `node`, or `None` if `node` was already
    /// removed.
    pub fn weak_in_neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let weak_in_neighbors = in_set.difference(self.out_list[node].as_ref().expect("`node` exists")).copied().collect::<FxHashSet<_>>();
            Some(weak_in_neighbors)
        } else {
            None
        }
    }

    /// Returns the set of weak outgoing neighbors of `node`, or `None` if `node` was already
    /// removed.
    pub fn weak_out_neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(out_set) = &self.out_list[node] {
            let weak_out_neighbors = out_set.difference(self.in_list[node].as_ref().expect("`node` exists")).copied().collect::<FxHashSet<_>>();
            Some(weak_out_neighbors)
        } else {
            None
        }
    }

    /// Returns a hashset of the strongly connected neighborhood of `node` that are inside `set`, 
    /// or None if `node` was already removed.
    pub fn strong_neighbors_in(&self, node: usize, set: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let out_set = &self.out_list[node].as_ref().expect("`node` exists");
            let mut strong_neighbors = in_set.intersection(out_set).copied().collect::<FxHashSet<_>>();
            strong_neighbors = strong_neighbors.intersection(set).copied().collect();
            Some(strong_neighbors)
        } else {
            None
        }
    }

    /// Returns the strongly connected closed neighborhood of `node` that are inside `set`, 
    /// or None if `node` was already removed.
    pub fn strong_closed_neighbors_in(&self, node: usize, set: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        if let Some(in_set) = &self.in_list[node] {
            let out_set = &self.out_list[node].as_ref().expect("`node` exists");
            let mut strong_neighbors = in_set.intersection(out_set).copied().collect::<FxHashSet<_>>();
            strong_neighbors.insert(node);
            strong_neighbors = strong_neighbors.intersection(set).copied().collect();
            Some(strong_neighbors)
        } else {
            None
        }
    }

    /// Returns a hashset of the out neighborhood of `node` that are also in `set` if `node` was 
    /// not already deleted.
    pub fn out_neighbors_in(&self, node: usize, set: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        self.out_list[node].as_ref().map(|outs| set.intersection(outs).copied().collect())
    }

    /// Returns a hashset of the in neighborhood of `node` that are also in `set` if `node` was
    /// not already deleted.
    pub fn in_neighbors_in(&self, node: usize, set: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        self.in_list[node].as_ref().map(|ins| set.intersection(ins).copied().collect())
    }

    /// Checks if `node` is a link node (has exectly two strong neighbors and no weak neighbors)
    /// and returns the neighbors if it is a link node and `None` otherwise.
    pub fn is_link_node(&self, node: usize) -> Option<[usize;2]> {
        if let Some(in_set) = &self.in_list[node] {
            if in_set.len() == 2 {
                let out_set = &self.out_list[node].as_ref().expect("`node` exists");
                if &in_set == out_set {
                    let mut sup = in_set.iter();
                    return Some([*sup.next().expect("has two elements"), *sup.next().expect("has two elements")]);
                }
            }
        }
        None 
    }

    /// Checks if `node` is a trip node (has exectly three strong neighbors and no weak neighbors)
    /// and returns the sorted neighbors if it is a trip node and `None` otherwise.
    pub fn is_trip_node(&self, node: usize) -> Option<[usize;3]> {
        if let Some(in_set) = &self.in_list[node] {
            if in_set.len() == 3 {
                let out_set = &self.out_list[node].as_ref().expect("`node` exists");
                if &in_set == out_set {
                    let mut sup = in_set.iter();
                    let mut arr = [*sup.next().expect("has three elements"), *sup.next().expect("has three elements"), *sup.next().expect("has three elements")];
                    arr.sort_unstable();
                    return Some(arr);
                }
            }
        }
        None 
    }

    /// Returns the first found link node together with its strong neighbors, or `None` if no link
    /// node exists.
    pub fn get_any_link_node_and_neighbors(&self) -> Option<(usize, [usize;2])> {
        for node in self.nodes() {
            if let Some(neighs) = self.is_link_node(node) {
                return Some((node, neighs))
            }
        }
        None
    }

    /// Returns the smallest set between the outgoing and incoming neighbors of `node` where neighbors in
    /// `set` are ignored, or None if `node` does not exist. 
    ///
    /// # Panics
    /// Panics if `node` is out of bounds.
    pub fn min_direct_neighbors_outside(&self, node: usize, set: &FxHashSet<usize>) -> Option<FxHashSet<usize>> {
        if let Some(ins) = &self.in_list[node] {
            let in_nei_wo: FxHashSet<usize> = ins.difference(set).copied().collect();
            let outs = &self.out_list[node].as_ref().expect("`node` exists");
            let out_nei_wo: FxHashSet<usize> = outs.difference(set).copied().collect();
            return if in_nei_wo.len() < out_nei_wo.len() {
                Some(in_nei_wo)
            } else {
                Some(out_nei_wo)
            }
        }
        None
    }

    /// Returns the minimum between the incoming neighbors of `node` and the outgoing neighbors of `node` if `node` was not already deleted.
    pub fn min_direct_neighbors(&self, node: usize) -> Option<FxHashSet<usize>> {
        if let Some(ins) = &self.in_list[node] {
            let outs = self.out_list[node].as_ref().expect("`node` exists");
            return if ins.len() < outs.len() {
                Some(ins.clone())
            } else {
                Some(outs.clone())
            }
        } 
        None
    }

    /// Returns the union of all out neighbors of `set`.
    /// Does not exclude neighbors in `set`.
    pub fn open_out_neighbors_of_set(&self, set: &FxHashSet<usize>) -> FxHashSet<usize> {
        set.iter().filter_map(move |node| self.out_list[*node].as_ref()).flatten().copied().collect()
    }

    /// Returns the amount of Nodes including the removed once.
    pub fn num_reserved_nodes(&self) -> usize {
        self.out_list.len()
    }

    /// Returns an iterator over all undeleted nodes.
    pub fn nodes(&self) -> impl Iterator<Item=usize> + '_ {
       self.out_list.iter()
           .enumerate()
           .filter_map(|(index, node)| {
               if node.is_some() {
                   Some(index)
               } else {
                   None
               }
           })
    }

    /// Returns the number of nodes in the graph.
    pub fn num_nodes(&self) -> usize {
        self.nodes().count()
    }

    /// Returns an iterator over all nodes that are not connected to any weak edges.
    pub fn only_strong_nodes(&self) -> impl Iterator<Item=usize> + '_ {
        self.nodes().filter(move |node| self.in_list[*node] == self.out_list[*node])
    }

    /// Checks if `node` exists.
    pub fn has_node(&self, node: usize) -> bool {
        self.out_list.len() > node && self.out_list[node].is_some() && self.in_list[node].is_some()
    }

    /// Returns the number of edges in the graph.
    pub fn num_edges(&self) -> usize {
        self.edges().count()
    }

    /// Checks if `edge` exists.
    pub fn has_edge(&self, edge: (usize, usize)) -> bool {
        self.out_list.len() > edge.0 && self.out_list[edge.0].as_ref().filter(|outs| outs.contains(&edge.1)).is_some()
    }

    /// Checks if `edge` exists and does not have a reciprocal twin.
    pub fn is_weak_edge(&self, edge: (usize, usize)) -> bool {
        self.out_list.len() > edge.0 && self.out_list[edge.0].as_ref().filter(|outs| outs.contains(&edge.1)).is_some() && self.in_list[edge.0].as_ref().filter(|ins| ins.contains(&edge.1)).is_none()
    }

    /// Checks if `node_a` is strongly connected to `node_b`.
    pub fn strongly_connected(&self, node_a: usize, node_b: usize) -> bool {
        if let Some(outs) = &self.out_list[node_a] {
            let ins = self.in_list[node_a].as_ref().expect("`node_a` exists");
            if outs.contains(&node_b) && ins.contains(&node_b) {
                return true
            }
        }
        false
    }

    /// Returns an iterator over all strong edges.
    pub fn strong_edges(&self) -> impl Iterator<Item=(usize, usize)> + '_ {
        self.in_list.iter()
            .enumerate()
            .filter(|(_,neighs)| neighs.is_some())
            .flat_map(move |(index, neighs)| {
                neighs.as_ref()
                    .expect("Due to filter")
                    .intersection(self.out_list[index].as_ref().expect("Due to filter"))
                    .copied()
                    .map(|nn| (index, nn))
                    .collect::<Vec<(usize, usize)>>()
            })
    }

    /// Returns an iterator over all edges.
    pub fn edges(&self) -> impl Iterator<Item=(usize, usize)> + '_ {
       self.out_list.iter()
           .enumerate()
           .filter(|(_,neighs)| neighs.is_some())
           .flat_map(|(index, neighs)| {
               neighs.as_ref()
                   .expect("Due to filter")
                   .iter()
                   .copied()
                   .map(|nn| (index, nn))
                   .collect::<Vec<(usize, usize)>>()
           })
    }

    /// Returns a `Vec` of edges from `from` to `to`.
    pub fn edges_from_to(&self, from: &FxHashSet<usize>, to: &FxHashSet<usize>) -> Vec<(usize, usize)> {
        from.iter()
            .copied()
            .filter(|node| self.out_list[*node].is_some())
            .flat_map(|node| {
                self.out_list[node].as_ref()
                    .expect("Due to filter()")
                    .iter()
                    .filter_map(move |neigbor| {
                        if to.contains(neigbor) {
                            Some((node, *neigbor))
                        } else {
                            None
                        }
                    })
            }).collect()
    }

    /// Checks if the nodes in `set` share a strong edge.
    pub fn has_strong_edge(&self, set: &FxHashSet<usize>) -> bool {
        let mut mut_set = set.clone();
        while mut_set.len() > 1 {
            let next = *mut_set.iter().next().expect("`mut_set` is not empty");
            mut_set.remove(&next);
            if let Some(strong_neighbors) = self.strong_neighbors(next) {
                if mut_set.intersection(&strong_neighbors).count() > 0 {
                    return true
                }
            }
        }
        false
    }

    /// Returns a `Vec` of weak edges from `from` to `to`.
    pub fn weak_edges_from_to(&self, from: &FxHashSet<usize>, to: &FxHashSet<usize>) -> Vec<(usize, usize)> {
        from.iter()
            .copied()
            .filter(|node| self.out_list[*node].is_some())
            .flat_map(|node| {
                self.out_list[node].as_ref()
                    .expect("Due to filter()")
                    .iter()
                    .filter_map(move |neigbor| {
                        if to.contains(neigbor) && 
                            !self.out_list[*neigbor].as_ref().expect("`neighbor` is neighbor of `node`").contains(&node) {
                            Some((node, *neigbor))
                        } else {
                            None
                        }
                    })
            }).collect()
    }

    /// Returns an iterator over all weak edges.
    pub fn weak_edges(&self) -> impl Iterator<Item=(usize, usize)> + '_ {
       self.out_list.iter()
           .enumerate()
           .filter_map(move |(node, neighs)| {
               if neighs.is_some() {
                    Some((node, neighs.as_ref().expect("is some").difference(self.in_list[node].as_ref().expect("`node` exists")).copied()))
               } else {
                   None
               }
           })
           .flat_map(|(index, neighs)| {
               neighs.map(|nn| (index, nn))
                   .collect::<Vec<(usize, usize)>>()
           })
    }

    /// Returns all edges between `set_a` and `set_b`.
    ///
    /// TODO: This implementation sucks!
    pub fn edges_between<'a: 'b, 'b>(&'a self, set_a: &'b FxHashSet<usize>, set_b: &'b FxHashSet<usize>) 
    -> impl Iterator<Item=(usize, usize)> + 'b {
        self.edges()
            .filter(move |(src, trg)| (set_a.contains(src) && set_b.contains(trg)) || 
                     (set_b.contains(src) && set_a.contains(trg)))
    }

    /// Returns all weak edges between `set_a` and `set_b`.
    ///
    /// TODO: This implementation sucks!
    pub fn weak_edges_between<'a: 'b, 'b>(&'a self, set_a: &'b FxHashSet<usize>, set_b: &'b FxHashSet<usize>) 
    -> impl Iterator<Item=(usize, usize)> + 'b {
        self.weak_edges()
            .filter(move |(src, trg)| (set_a.contains(src) && set_b.contains(trg)) || 
                     (set_b.contains(src) && set_a.contains(trg)))
    }

    /// Returns the out degree of `node` or None if `node` was deleted. 
    ///
    /// # Panics
    /// Panics if node id is out of bounds.
    pub fn out_degree(&self, node: usize) -> Option<usize> {
        self.out_list[node].as_ref().map(|outs| outs.len())
    }

    /// Returns the in degree of `node` or None if `node` was deleted. 
    ///
    /// # Panics
    /// Panics if node id is out of bounds.
    pub fn in_degree(&self, node: usize) -> Option<usize> {
        self.in_list[node].as_ref().map(|ins| ins.len())
    }

    /// Returns the minimal weak degree of `node` or `None` if `node` does not exist.
    pub fn min_weak_degree(&self, node: usize) -> Option<usize> {
        if let Some(ins) = &self.in_list[node] {
            let outs = self.out_list[node].as_ref().expect("`node` does not exist");
            let weak_out_deg = outs.difference(ins).count();
            let weak_in_deg = ins.difference(outs).count();
            Some(min(weak_in_deg, weak_out_deg))
        } else {
            None
        }
    }

    /// Returns the degree of `node` or None if `node` was deleted. 
    ///
    /// # Panics
    /// Panics if node id is out of bounds.
    pub fn degree(&self, node: usize) -> Option<usize> {
        if let Some(in_degree) = self.in_degree(node) {
            self.out_degree(node).map(|out_degree| in_degree + out_degree)
        } else {
            None
        }
    }

    /// Returns the strong degree (amount of nodes that that are connected via both an incoming-
    /// and an outgoing edge) of `node` or None if `node` was deleted. 
    ///
    /// # Panics 
    /// Panics if `node` is out of bounds.
    pub fn strong_degree(&self, node: usize) -> Option<usize> {
        self.strong_neighbors(node).as_ref().map(|neighs| neighs.len())
    }

    /// Returns the minimum between outgoing and incoming degree of `node` or None if `node` was
    /// deleted.
    ///
    /// # Panics
    /// Panics if node id is out of bounds.
    pub fn min_direct_degree(&self, node: usize) -> Option<usize> {
        if let Some(in_degree) = self.in_degree(node) {
            self.out_degree(node).map(|out_degree| min(in_degree,out_degree))
        } else {
            None
        }
    }


    /// Returns the minimum between the outgoing and incoming degree of `node` where neighbors in
    /// `set` are ignored, or None if `node` does not exist. 
    ///
    /// # Panics
    /// Panics if `node` is out of bounds.
    pub fn min_direct_degree_outside(&self, node: usize, set: &FxHashSet<usize>) -> Option<usize> {
        if let Some(ins) = &self.in_list[node] {
            let in_deg_wo = ins.difference(set).count();
            let outs = &self.out_list[node].as_ref().expect("`node` exists");
            let out_deg_wo = outs.difference(set).count();
            return Some(min(out_deg_wo, in_deg_wo))
        }
        None
    }

    /// Weight function proposed by Cai et al. "A search algorithm for computing minimum feedback
    /// vertex set of a directed graph."
    ///
    /// Needs only a single weight, second can be 0.
    pub fn cai_weight(&self, node: usize, weights: (f64, f64)) -> Option<f64> {
        if let Some(incoming) = self.in_degree(node) {
            self.out_degree(node).map(|outgoing| incoming as f64 + outgoing as f64 - weights.0 * (incoming as f64 - outgoing as f64).abs())
        } else {
            None
        }
    }

    /// Weight function proposed by Lin et al. "On computing the minimum feedback vertex set of a 
    /// directed graph by contraction operations"
    pub fn lin_weight(&self, node: usize, weights: (f64, f64)) -> Option<f64> {
        if let Some(strong) = self.strong_degree(node) {
            self.degree(node).map(|degree| weights.0 * strong as f64 + weights.1 * degree as f64)
        } else {
            None
        }
    }

    /// Returns the node with the highest weight computed by `fun`, or None if no node exists. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_max_weight_node(&self, fun: &dyn Fn(&Self, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64)) -> Option<usize> {
        let mut cur_max_weight = 0f64;
        let mut cur_max_node: Option<usize> = None;
        for node in self.nodes() {
            let weight = fun(self, node, weights).expect("This node exists, since it is in .nodes()");
            if weight > cur_max_weight {
                cur_max_weight = weight;
                cur_max_node = Some(node);
            }
        }
        cur_max_node
    }

    /// Returns the node with the lowest weight computed by `fun`, or None if no node exists. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_min_weight_node(&self, fun: &dyn Fn(&Self, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64)) -> Option<usize> {
        let mut cur_min_weight = None;
        let mut cur_min_node: Option<usize> = None;
        for node in self.nodes() {
            let weight = fun(self, node, weights).expect("This node exists, since it is in .nodes()");
            if let Some(min_weight) = cur_min_weight {
                if weight < min_weight {
                    cur_min_weight = Some(weight);
                    cur_min_node = Some(node);
                }
            } else {
                cur_min_node = Some(node);
                cur_min_weight = Some(weight);
            }
        }
        cur_min_node
    }

    /// Returns the node with the lowest weight computed by `fun`, or None if no node is in `set`
    /// or `set` is not completely in `self`. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_min_weight_node_within(&self, fun: &dyn Fn(&Self, usize, (f64, f64)) -> Option<f64>, weights: (f64, f64), set: &FxHashSet<usize>) -> Option<usize> {
        let mut cur_min_weight = None;
        let mut cur_min_node: Option<usize> = None;
        for node in set {
            if let Some(weight) = fun(self, *node, weights){
                if let Some(min_weight) = cur_min_weight {
                    if weight < min_weight {
                        cur_min_weight = Some(weight);
                        cur_min_node = Some(*node);
                    }
                } else {
                    cur_min_node = Some(*node);
                    cur_min_weight = Some(weight);
                }
            } else {
                return None
            }
        }
        cur_min_node
    }

    /// Returns the node with the highest degree, or None if no node exists. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_max_degree_node(&self) -> Option<usize> {
        let mut cur_max_degree: usize = 0;
        let mut cur_max_node: Option<usize> = None;
        for node in self.nodes() {
            let degree = self.degree(node).expect("This node exists, since it is in .nodes()");
            if degree > cur_max_degree {
                cur_max_degree = degree;
                cur_max_node = Some(node);
            }
        }
        cur_max_node
    }

    //pub fn find_reducible_double_cycles(&self) -> HashSet<Vec<usize>> {
    //    let mut used = FxHashSet::default();
    //    let mut visited = FxHashSet::default();
    //    // Let nodes with st_deg == 2 and wk_deg == 0 be called link nodes
    //    for node in self.nodes() {
    //        let mut queue = Vec::new();
    //        if used.contains(&node) {
    //            continue;
    //        }
    //        used.insert(node);
    //        queue.extend(self.strong_neighbors(node).iter());
    //        while !queue.is_empty() {
    //            let 
    //            }
    //    }
    //    return HashSet::new()
    //    // For start in nodes:
    //    //  If alread marked as visited, skip 
    //    //  Mark start as visited
    //    //  put all strong neighbors in queue 
    //    //  for neigh in queue: 
    //    //      if neigh.is_visited: 
    //    //          backtrace and check if we found a good cycle 
    //    //      mark neigh as visited.
    //    //      add all strong neighbors of neigh into queue 
    //    //
    //    // If we go back we have to unmark visited nodes!!!
    //    //      
    //}

    /// Returns the node with the highest minimum direct degree, or None if no node with at least one incoming- and one outgoing neighbor exists. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_max_min_direct_degree_node(&self) -> Option<usize> {
        let mut cur_max_degree: usize = 0;
        let mut cur_max_node: Option<usize> = None;
        for node in self.nodes() {
            let degree = self.min_direct_degree(node).expect("This node exists, since it is in .nodes()");
            if degree > cur_max_degree {
                cur_max_degree = degree;
                cur_max_node = Some(node);
            }
        }
        cur_max_node
    }

    /// Returns the node with the lowest minimum direct degree and the respective neighbors, or None if no such node exists. 
    ///
    /// # Panics
    /// Panics if `self` is broken.
    pub fn get_min_min_direct_degree_node_and_neighbors(&self) -> Option<(usize, FxHashSet<usize>)> {
        let mut cur_min_degree: Option<usize> = None;
        let mut cur_min_pair: Option<(usize, FxHashSet<usize>)> = None;
        for node in self.nodes() {
            let min_neighbors = self.min_direct_neighbors(node).expect("This node exists, since it is in .nodes()");
            if let Some(min_degree) = cur_min_degree {
                if min_neighbors.len() < min_degree {
                    cur_min_degree = Some(min_neighbors.len());
                    cur_min_pair = Some((node, min_neighbors));
                }
            } else {
                cur_min_degree = Some(min_neighbors.len());
                cur_min_pair = Some((node, min_neighbors));
            }
        }
        cur_min_pair
    }

    /// Returns the node with the highest strong degree (see .strong_degree) and its strong
    /// neighborhood, or None if no node with strong neighbors exist.
    pub fn get_max_strong_degree_node_and_neighbors(&self) -> Option<(usize, FxHashSet<usize>)> {
        let mut cur_max_degree: usize = 0; // To save some time?
        let mut cur_max_set: Option<(usize, FxHashSet<usize>)> = None;
        for node in self.nodes() {
            let strong_neighborhood = self.strong_neighbors(node).expect("`node` exists, since it is in .nodes()");
            let degree = strong_neighborhood.len();
            if degree > cur_max_degree {
                // TODO: less cloning to save more time?
                cur_max_set = Some((node, strong_neighborhood.clone()));
                cur_max_degree = degree;
            }
        }
        cur_max_set
    }

    /// Checks if `set` is a strongly connected cluster.
    pub fn is_cluster(&self, set: &FxHashSet<usize>) -> bool {
        for node in set {
            if self.in_list[*node].is_none() {
                return false
            }
            if !set.is_subset(&self.closed_strong_neighbors(*node).expect("`node` exists")) {
                return false
            }
        }
        true
    }

    /// Heuristic to check if it makes sense to treat graph as a vc instance.
    pub fn is_vcable(&self) -> bool {
        if self.num_nodes() < 50 {
            return false
        }
        if ((self.strong_edges().count() * 2) as f64) < self.num_edges() as f64 * 0.80 {
            return false
        }
        return true
    }

    /// Checks if the graph holds the edge `(src, trg)`.
    pub fn edge_exists(&self, (src, trg): &(usize, usize)) -> bool {
        if let Some(src_out_list) = &self.out_list[*src] {
            return src_out_list.contains(trg)
        }
        false
    }

    /// Returns a maximal matching of edges within `set_a` or between `set_a` and `set_b` if
    /// `set_b` is given and the set of unmatched nodes within `set_a`.
    ///
    /// TODO: Currently unknown behaviour if any nodes in `set_a` or in `set_b` is not a strong
    /// node.
    ///
    /// TODO: Needs testing
    pub fn strong_max_matching_between(&self, set_a: &FxHashSet<usize>, set_b: &Option<FxHashSet<usize>>) -> (HashSet<(usize, usize)>, FxHashSet<usize>) {
        let mut matched_edges: HashSet<(usize, usize)> = HashSet::new();
        let mut unmatched_nodes = FxHashSet::default();
        if let Some(mut ex_set_b) = set_b.clone() {
            let mut ex_set_a = set_a.clone();
            while !(ex_set_a.is_empty() || ex_set_b.is_empty()) {
                let ele = ex_set_a.iter().next().cloned().expect("`ex_set_a` is not empty");
                ex_set_a.remove(&ele);
                // all nodes are strong nodes:ele
                if let Some(&neigh) =  self.in_list[ele].as_ref().expect("`ele` should exist").intersection(&ex_set_b).next().clone() {
                    matched_edges.insert((ele, neigh));
                    ex_set_b.remove(&neigh);
                } else {
                    unmatched_nodes.insert(ele);    
                }
            }
            unmatched_nodes.extend(ex_set_a.iter());
        } else {
            let mut ex_set_a = set_a.clone();
            while !ex_set_a.is_empty() {
                let ele = ex_set_a.iter().next().cloned().expect("`ex_set_a` is not empty");
                ex_set_a.remove(&ele);
                // all nodes are strong nodes:ele
                if let Some(&neigh) =  self.in_list[ele].as_ref().expect("`ele` should exist").intersection(&ex_set_a).next() {
                    matched_edges.insert((ele, neigh));
                    ex_set_a.remove(&neigh);
                } else {
                    unmatched_nodes.insert(ele);    
                }
            }
        }
        (matched_edges, unmatched_nodes)
    }

    /// Find all nodes reachable from `node`. 
    ///
    /// # Panics
    /// Can panic if invalid nodes are in `out_list`, or `node` is an invalid id or deleted.
    pub fn undirected_reachable(&self, node: usize) -> FxHashSet<usize> {
        let mut reached: FxHashSet<usize> = FxHashSet::default();
        let mut queue: Vec<usize> = vec![node];
        while !queue.is_empty() {
            let current = queue.pop().expect("`queue` is not empty.");
            if reached.contains(&current) {
                continue;
            }
            reached.insert(current);
            queue.append(&mut self.neighbors(current).expect("This node has to exist").into_iter().collect::<Vec<usize>>());
        }
        reached
    }

    /// Finds all nodes reachable from `node` by directed edges. 
    ///
    /// # Panics
    /// Can panic if invalid nodes are in `out_list`, or `node` is an invalid id or deleted.
    pub fn directed_reachable(&self, node: usize) -> FxHashSet<usize> {
        let mut reached: FxHashSet<usize> = FxHashSet::default();
        let mut queue: Vec<usize> = vec![node];
        while !queue.is_empty() {
            let current = queue.pop().expect("`queue` is not empty.");
            if reached.contains(&current) {
                continue;
            }
            reached.insert(current);
            let out_neighbors = self.out_neighbors(current).clone().expect("This node has to exist");
            queue.append(&mut out_neighbors.into_iter().collect::<Vec<usize>>());
        }
        reached
    }

    /// Checks if any node in `trg` can be reached by weak directed edges from any node in `src`.
    ///
    /// # Panics 
    /// Panics if any node in `src` does not exist.
    pub fn weak_path_exists_between(&self, src: &FxHashSet<usize>, trg: &FxHashSet<usize>)
        -> bool {
        let mut queue: Vec<usize> = src.iter().copied().collect();
        let mut marked: FxHashSet<usize> = FxHashSet::default();
        while !queue.is_empty() {
            let current = queue.pop().expect("`queue` is not empty.");
            if trg.contains(&current) {
                return true
            }
            if marked.contains(&current) {
                continue;
            }
            marked.insert(current);
            let weak_out_neighbors = self.weak_out_neighbors(current).clone().expect("`current` should exist");
            queue.append(&mut weak_out_neighbors.into_iter().collect::<Vec<usize>>());
        }
        false
    }

    /// Checks if a graph `self`, that has only strongly connected components, is disconnected.
    ///
    /// Does not work on other graphs.
    pub fn scc_graph_is_disconnected(&self) -> bool {
        let num_nodes = self.num_nodes();
        if num_nodes == 0 {
            return false
        }
        let reachable_for_some = self.directed_reachable(self.nodes().next().expect("`self` is not empty")); 
        reachable_for_some.len() != num_nodes
    }

    /// Counts the node disjunct cycles starting at `node` and return them.
    pub fn count_petals_give_petals(&self, node: usize) -> (usize, FxHashSet<usize>) {
        let daisy_petals = self.strong_neighbors(node).expect("`node` exists");
        let mut srcs = self.out_neighbors(node).as_ref().expect("`node` exitsts").clone();
        srcs = srcs.difference(&daisy_petals).copied().collect();
        let mut trgs = self.in_neighbors(node).as_ref().expect("`node` exitsts").clone();
        trgs = trgs.difference(&daisy_petals).copied().collect();
        let mut num_petals = daisy_petals.len();
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for src in &srcs {
            // All sources in `srcs` besides `src` and all nodes in `daisy_petals` can not be used
            // for any other petal.
            let mut preds: FxHashMap<usize, usize> = srcs.iter().copied()
                .filter_map(|s| if s==*src {None} else{Some((s, s))}).collect();
            preds.extend(daisy_petals.iter().map(|d| (d,d)));
            let mut queue = VecDeque::new();
            queue.push_back((*src, *src));
            while !queue.is_empty() {
                let (next, pref) = queue.pop_front().expect("`queue` not empty");
                if preds.contains_key(&next) {
                    // If `next` is marked.
                    continue;
                }
                if taken.contains_key(&next) {
                    if !taken.contains_key(&pref) {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        preds.insert(next, pref);
                        continue;
                    } else {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        self.out_neighbors(next).as_ref().expect("`src` exists").iter().filter(|n| **n != pref).for_each(|n| queue.push_back((*n, next)));
                        preds.insert(next, pref);
                        continue;
                    }
                } 
                if trgs.contains(&next) {
                    num_petals += 1;
                    // backtrack
                    let mut cur = next;
                    let mut pre = pref;
                    taken.insert(cur, pre);
                    let mut before = None;
                    while pre != *src {
                        cur = pre;
                        pre = *preds.get(&cur).expect("`cur` should have been traversed");
                        if taken.contains_key(&cur) {
                            if before.is_some() {
                                taken.remove(&cur);
                            }
                            before = Some((cur,pre));
                        } else {
                            if let Some((p1, p2)) = before {
                                assert_ne!(p1, p2);
                                taken.insert(p1, p2);
                                before = None;
                            }
                            taken.insert(cur, pre);
                        }
                    }
                    // `taken` can be used to find the left over graph.
                    taken.insert(*src, *src);
                    break
                }
                self.out_neighbors(next).as_ref().expect("`src` exists").iter().for_each(|n| queue.push_back((*n, next)));
                preds.insert(next, pref);
            }
        }
        // Left over graph can be build by removing all nodes in 
        // `taken`, `daisy_petals` and `node`
        let mut petal_nodes: FxHashSet<_> = taken.keys().copied().chain(daisy_petals.iter().copied()).collect();
        petal_nodes.insert(node);
        (num_petals, petal_nodes)
    }

    /// Counts the node disjunct cycles starting at `node` and creates a clone of `self` with all
    /// nodes of the cycles removed.
    pub fn count_petals_left_over(&self, node: usize) -> (usize, Self) {
        let daisy_petals = self.strong_neighbors(node).expect("`node` exists");
        let mut srcs = self.out_neighbors(node).as_ref().expect("`node` exitsts").clone();
        srcs = srcs.difference(&daisy_petals).copied().collect();
        let mut trgs = self.in_neighbors(node).as_ref().expect("`node` exitsts").clone();
        trgs = trgs.difference(&daisy_petals).copied().collect();
        let mut num_petals = daisy_petals.len();
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for src in &srcs {
            // All sources in `srcs` besides `src` and all nodes in `daisy_petals` can not be used
            // for any other petal.
            let mut preds: FxHashMap<usize, usize> = srcs.iter().copied()
                .filter_map(|s| if s==*src {None} else{Some((s, s))}).collect();
            preds.extend(daisy_petals.iter().map(|d| (d,d)));
            let mut queue = VecDeque::new();
            queue.push_back((*src, *src));
            while !queue.is_empty() {
                let (next, pref) = queue.pop_front().expect("`queue` not empty");
                if preds.contains_key(&next) {
                    // If `next` is marked.
                    continue;
                }
                if taken.contains_key(&next) {
                    if !taken.contains_key(&pref) {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        preds.insert(next, pref);
                        continue;
                    } else {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        self.out_neighbors(next).as_ref().expect("`src` exists").iter().filter(|n| **n != pref).for_each(|n| queue.push_back((*n, next)));
                        preds.insert(next, pref);
                        continue;
                    }
                } 
                if trgs.contains(&next) {
                    num_petals += 1;
                    // backtrack
                    eprintln!("Found petal");
                    let mut cur = next;
                    let mut pre = pref;
                    eprintln!("Edge: {} {}", cur, pre);
                    taken.insert(cur, pre);
                    let mut before = None;
                    while pre != *src {
                        cur = pre;
                        pre = *preds.get(&cur).expect("`cur` should have been traversed");
                        if taken.contains_key(&cur) {
                            if before.is_some() {
                                taken.remove(&cur);
                            }
                            before = Some((cur,pre));
                        } else {
                            if let Some((p1, p2)) = before {
                                assert_ne!(p1, p2);
                                taken.insert(p1, p2);
                                before = None;
                            }
                            eprintln!("Edge: {} {}", cur, pre);
                            taken.insert(cur, pre);
                        }
                    }
                    // `taken` can be used to find the left over graph.
                    taken.insert(*src, *src);
                    break
                }
                self.out_neighbors(next).as_ref().expect("`src` exists").iter().for_each(|n| queue.push_back((*n, next)));
                preds.insert(next, pref);
            }
        }
        // Left over graph can be build by removing all nodes in 
        // `taken`, `daisy_petals` and `node`
        let mut clone = self.clone();
        clone.remove_nodes(taken.keys().copied());
        clone.remove_nodes(daisy_petals.iter().copied());
        clone.remove_node(node);
        (num_petals, clone)
    }

    /// Uses max-flow min-cut to compute the maximal amount of node disjunct cycles starting at
    /// `node`.
    pub fn count_petals(&self, node: usize) -> usize {
        let daisy_petals = self.strong_neighbors(node).expect("`node` exists");
        let mut srcs = self.out_neighbors(node).as_ref().expect("`node` exitsts").clone();
        srcs = srcs.difference(&daisy_petals).copied().collect();
        let mut trgs = self.in_neighbors(node).as_ref().expect("`node` exitsts").clone();
        trgs = trgs.difference(&daisy_petals).copied().collect();
        let mut num_petals = daisy_petals.len();
        let mut taken: HashMap<usize, usize> = HashMap::new();
        for src in &srcs {
            // All sources in `srcs` besides `src` and all nodes in `daisy_petals` can not be used
            // for any other petal.
            let mut preds: FxHashMap<usize, usize> = srcs.iter().copied()
                .filter_map(|s| if s==*src {None} else{Some((s, s))}).collect();
            preds.extend(daisy_petals.iter().map(|d| (d,d)));
            let mut queue = VecDeque::new();
            queue.push_back((*src, *src));
            while !queue.is_empty() {
                let (next, pref) = queue.pop_front().expect("`queue` not empty");
                if preds.contains_key(&next) {
                    // If `next` is marked.
                    continue;
                }
                if taken.contains_key(&next) {
                    if !taken.contains_key(&pref) {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        preds.insert(next, pref);
                        continue;
                    } else {
                        queue.push_back((*taken.get(&next).expect("`next` exists"), next));
                        // keep in mind reverse if this route is taken
                        self.out_neighbors(next).as_ref().expect("`src` exists").iter().filter(|n| **n != pref).for_each(|n| queue.push_back((*n, next)));
                        preds.insert(next, pref);
                        continue;
                    }
                } 
                if trgs.contains(&next) {
                    num_petals += 1;
                    // backtrack
                    let mut cur = next;
                    let mut pre = pref;
                    taken.insert(cur, pre);
                    let mut before = None;
                    while pre != *src {
                        cur = pre;
                        pre = *preds.get(&cur).expect("`cur` should have been traversed");
                        if taken.contains_key(&cur) {
                            if before.is_some() {
                                taken.remove(&cur);
                            }
                            before = Some((cur,pre));
                        } else {
                            if let Some((p1, p2)) = before {
                                assert_ne!(p1, p2);
                                taken.insert(p1, p2);
                                before = None;
                            }
                            taken.insert(cur, pre);
                        }
                    }
                    // `taken` can be used to find the left over graph.
                    taken.insert(*src, *src);
                    break
                }
                self.out_neighbors(next).as_ref().expect("`src` exists").iter().for_each(|n| queue.push_back((*n, next)));
                preds.insert(next, pref);
            }
        }
        // Left over graph can be build by removing all nodes in 
        // `taken` and `daisy_petals`
        num_petals
    }

    /// Finds a set of n nodes, such that their removal would split the graph into at least n+1
    /// strongly connected components.
    pub fn find_split_set(&self) -> FxHashSet<usize> {
        // get all sccs 
        let subs = self.split_into_connected_components_alt();
        let mut cut_vertices: FxHashSet<usize> = FxHashSet::default();
        for mut sub in subs {
            if let Some(cut_vs) = sub.split_scc_recursively() {
                cut_vertices.extend(cut_vs.into_iter());
            }
        }
        return cut_vertices;
    }

}

// Dynamic operations.
impl Digraph {

    /// Looks for a node that splits `self` into two or more strongly connected components and
    /// repeats recursively until no more node can be found. Returns the set of cut nodes, or
    /// `None` if none could have been found.
    fn split_scc_recursively(&mut self) -> Option<FxHashSet<usize>> {
        let mut nodes: Vec<usize> = self.nodes().collect();
        nodes.sort_unstable_by_key(|n| self.degree(*n));
        while !nodes.is_empty() {
            let tn = nodes.pop().expect("`nodes` is not empty");
            let mut clone = self.clone();
            clone.remove_node(tn);
            if let Some(sccs) = clone.reduce_to_sccs() {
                let mut cut_vertices: FxHashSet<usize> = vec![tn].into_iter().collect();
                if sccs.len() > 1 {
                    let subs = clone.split_into_connected_components_alt();
                    for mut sub in subs {
                        if let Some(cut_vs) = sub.split_scc_recursively(){
                            cut_vertices.extend(cut_vs.into_iter());
                        }
                    }
                    return Some(cut_vertices);
                }
            }
        }
        None
    }

    /// Removes the node `node` and all adjacent edges from the graph.
    /// Returns the a tuple of incoming and outgoing neighbors if the node was removed
    /// successfully. Returns None if the node did not exist.
    ///
    /// # Panics
    /// Panics if the node index is out of bounds or the graph is broken.
    pub fn remove_node(&mut self, node: usize) -> Option<(FxHashSet<usize>, FxHashSet<usize>)> {
        if let Some(in_neighbors) = self.in_list[node].take() {
            for in_neighbor in &in_neighbors {
                if *in_neighbor != node {
                    self.out_list[*in_neighbor].as_mut().unwrap().remove(&node);
                }
            }
            if let Some(out_neighbors) = self.out_list[node].take() {
                for out_neighbor in &out_neighbors {
                    if *out_neighbor != node {
                        self.in_list[*out_neighbor].as_mut().unwrap().remove(&node);
                    }
                }
                Some((in_neighbors, out_neighbors))
            } else {
                // This should never happen.
                panic!();
            }
        } else {
            None
        }
    }

    /// Removes `nodes` from the graph
    ///
    /// # Panics
    /// Panics if a node index is invalid.
    pub fn remove_nodes<I: IntoIterator<Item=usize>>(&mut self, nodes: I) {
        for node in nodes {
            self.remove_node(node);
        }
    }

    /// Reinsert a node that was deleted earlier.
    /// Returns false if the node was not deleted properly.
    pub fn reinsert_node(&mut self, node: usize, in_neighbors: FxHashSet<usize>,
                         out_neighbors: FxHashSet<usize>) -> bool {
        if self.in_list[node].is_none() {
            for in_neighbor in &in_neighbors {
                if *in_neighbor != node {
                    self.out_list[*in_neighbor].as_mut().expect("It was there as `node` was deleted").insert(node);
                }
            }
            self.in_list[node] = Some(in_neighbors);
        } else {
            return false
        }
        if self.out_list[node].is_none() {
            for out_neighbor in &out_neighbors {
                if *out_neighbor != node {
                    self.in_list[*out_neighbor].as_mut().expect("It was there as `node` was deleted").insert(node);
                }
            }
            self.out_list[node] = Some(out_neighbors);
        } else {
            return false // This means that the graph is broken.
        }
        true
    }

    /// Removes the directed edge `(src, trg)` from the graph.
    /// Returns if edge was removed
    ///
    /// # Panics
    /// Panics if a node index is invalid TODO: deleted also invalid?.
    pub fn remove_edge(&mut self, (src, trg): &(usize, usize)) -> bool {
        let e = self.out_list[*src].as_mut().unwrap().remove(trg);
        self.in_list[*trg].as_mut().unwrap().remove(src);
        e
    }

    /// Adds the directed edge `(src, trg)` to the graph if not already present.
    ///
    /// # Panics
    /// Panics if a node index is invalid TODO: deleted also invalid?.
    pub fn add_edge(&mut self, (src, trg): (usize, usize)) {
        self.out_list[src].as_mut().unwrap().insert(trg);
        self.in_list[trg].as_mut().unwrap().insert(src);
    }

    /// Adds the directed edge `(src, trg)` to the graph if not already present.
    /// Returns true if the edge was added.
    ///
    /// # Panics
    /// Panics if a node index is invalid 
    pub fn add_edge_checked(&mut self, (src, trg): (usize, usize)) -> bool {
        if let Some(out_list) = self.out_list[src].as_mut() {
            if let Some(in_list) = self.in_list[trg].as_mut() {
                if !out_list.contains(&trg) {
                    out_list.insert(trg);
                    in_list.insert(src);
                    return true
                }
            }
        }
        false
    }

    /// Removes directed edges of the form `(src, trg)` from the graph
    ///
    /// # Panics
    /// Panics if a node index is invalid.
    pub fn remove_edges<I: IntoIterator<Item=(usize, usize)>>(&mut self, edges: I) {
        for edge in edges {
            self.remove_edge(&edge);
        }
    }

    /// Adds directed edges of the form `(src, trg)` to the graph
    ///
    /// # Panics
    /// Panics if a node index is invalid.
    pub fn add_edges<I: IntoIterator<Item=(usize, usize)>>(&mut self, edges: I) {
        for edge in edges {
            self.add_edge(edge);
        }
    }

    /// Removes the directed edge `(src, trg)` and adds the directed edge `(trg, src)` to `self`
    /// Adds nothing if no edge was removed.
    ///
    /// # Panics
    /// Panics at invalid node index
    pub fn rev_edge(&mut self, (src, trg): (usize, usize)) {
        if self.remove_edge(&(src, trg)) {
            self.add_edge((trg, src));
        }
    }

    /// Reverses directed edges of the form `(src, trg)`.
    /// If an edge is not present, nothing happens.
    ///
    /// # Panics
    /// Panics if a node index is invalid.
    pub fn rev_edges<I: IntoIterator<Item=(usize, usize)>>(&mut self, edges: I) {
        for edge in edges {
            self.rev_edge(edge);
        }
    }


    /// Merges `node` with it's only incoming neighbor `into`.
    /// Returns the outgoing neighbors of `node` and the outgoing neighbors of `into` that are also
    /// outgoing neighbors of `node`.
    /// Returns None if `node` or `into` do not exist, or if `into` is not the only incoming neighbor
    /// of `node`.
    ///
    /// # Panics
    /// Panics if any index is out of bounds
    pub fn merge_into_back(&mut self, node: usize, into: usize) -> Option<(FxHashSet<usize>, FxHashSet<usize>)>{
        if self.in_degree(node).is_none() || self.in_degree(node).expect("First condition would be true") != 1 {
            return None
        }
        if !self.remove_edge(&(into, node)) {
            return None
        }
        let doubles: FxHashSet<usize> = self.out_neighbors(node)
            .as_ref()
            .expect("The existance of `node` was checked earlier")
            .clone()
            .iter()
            .copied()
            .filter(|neigh| {
                let in_list = self.in_list[*neigh].as_mut().expect("This node should exist");
                in_list.remove(&node);
                if in_list.contains(&into) {
                    true
                } else {
                    in_list.insert(into);
                    false
                }
            }).collect();
        let outs: FxHashSet<usize> = self.out_neighbors(node).as_ref().expect("The existance of `node` was checked earlier").clone();
        if let Some(into_outs) = self.out_list[into].clone()  {
            self.out_list[into] = Some(into_outs.union(&outs).copied().collect());
        } else {
            // TODO: Better to panic here, since the graph was already altered.
            return None
        }
        self.in_list[node] = None;
        self.out_list[node] = None;
        Some((outs, doubles))
    }

    /// Merges `node` with it's only outgoing neighbor `into`.
    /// Returns the incoming neighbors of `node` and the incoming neighbors of `into` that are also
    /// incoming neighbors of `node`.
    /// Returns None if `node` or `into` do not exist, or if `into` is not the only outgoing neighbor 
    /// of `node`.
    ///
    /// # Panics
    /// Panics if any index is out of bounds
    pub fn merge_into_front(&mut self, node: usize, into: usize) -> Option<(FxHashSet<usize>, FxHashSet<usize>)>{
        if self.out_degree(node).is_none() || self.out_degree(node).expect("First condition would be true") != 1 {
            return None
        }
        if !self.remove_edge(&(node, into)) {
            return None
        }
        let doubles: FxHashSet<usize> = self.in_neighbors(node)
            .as_ref()
            .expect("The existance of `node` was checked earlier")
            .clone()
            .iter()
            .copied()
            .filter(|neigh| {
                let out_list = self.out_list[*neigh].as_mut().expect("This node should exist");
                out_list.remove(&node);
                if out_list.contains(&into) {
                    true
                } else {
                    out_list.insert(into);
                    false
                }
            }).collect();
        let ins: FxHashSet<usize> = self.in_neighbors(node).as_ref().expect("The existance of `node` was checked earlier").clone();
        if let Some(into_ins) = self.in_list[into].clone()  {
            self.in_list[into] = Some(into_ins.union(&ins).copied().collect());
        } else {
            // TODO: Better to panic here, since the graph was already altered.
            return None
        }
        self.in_list[node] = None;
        self.out_list[node] = None;
        Some((ins, doubles))
    }

    /// Merges `node` into `into`, ignoring edges between those two nodes. 
    /// Returns a tuple with incoming and outgoing neighbors of only `node` and a tuple with incoming
    /// and outgoing neighbors of both `node` and `into`. Or `None` if something went wrong.
    ///
    /// # Panics
    /// Panics if any index is out of bounds
    pub fn complete_merge(&mut self, node: usize, into: usize)
        -> Option<((FxHashSet<usize>, FxHashSet<usize>), 
                   (FxHashSet<usize>, FxHashSet<usize>), 
                   HashSet<(usize, usize)>)> {
        let mut betweens = HashSet::new();
        if self.edge_exists(&(node, into)) {
            self.remove_edge(&(node, into));
            betweens.insert((node, into));
        } 
        if self.edge_exists(&(into, node)) {
            self.remove_edge(&(into, node));
            betweens.insert((into, node));
        } 
        // Removes `into` from it's neighbors, adds `into` to them, and records neighbors of both
        // `into` and `node`.
        let in_doubles: FxHashSet<usize> = self.in_neighbors(node)
            .as_ref()
            .expect("The existance of `node` was checked earlier")
            .clone()
            .iter()
            .copied()
            .filter(|neigh| {
                let out_list = self.out_list[*neigh].as_mut().expect("This node should exist");
                out_list.remove(&node);
                if out_list.contains(&into) {
                    true
                } else {
                    out_list.insert(into);
                    false
                }
            }).collect();
        let out_doubles: FxHashSet<usize> = self.out_neighbors(node)
            .as_ref()
            .expect("The existance of `node` was checked earlier")
            .clone()
            .iter()
            .copied()
            .filter(|neigh| {
                let in_list = self.in_list[*neigh].as_mut().expect("This node should exist");
                in_list.remove(&node);
                if in_list.contains(&into) {
                    true
                } else {
                    in_list.insert(into);
                    false
                }
            }).collect();
        let ins: FxHashSet<usize> = self.in_neighbors(node).as_ref().expect("The existance of `node` was checked earlier").clone();
        if let Some(into_ins) = self.in_list[into].clone()  {
            self.in_list[into] = Some(into_ins.union(&ins).copied().collect());
        } else {
            // TODO: Probably better to panic here, since the graph was already altered.
            return None
        }
        let outs: FxHashSet<usize> = self.out_neighbors(node).as_ref().expect("The existance of `node` was checked earlier").clone();
        if let Some(into_outs) = self.out_list[into].clone()  {
            self.out_list[into] = Some(into_outs.union(&outs).copied().collect());
        } else {
            // TODO: Probably better to panic here, since the graph was already altered.
            return None
        }
        self.in_list[node] = None;
        self.out_list[node] = None;
        Some(((ins, outs),(in_doubles, out_doubles), betweens))
    }

    /// Merges a list of nodes into another node. Ignores loops.
    /// Returns the lists of old incoming and old outgoing neighbors. 
    /// Throws an error if any of the given nodes are not in `self`.
    pub fn big_merge(&mut self, into: usize, from: Vec<usize>) -> Result<(Vec<FxHashSet<usize>>, Vec<FxHashSet<usize>>), ProcessingError> {
        let mut all_outs = Vec::new();
        let mut all_ins = Vec::new();
        if let Some(ins) = self.in_neighbors(into) {
            all_ins.push(ins.clone());
            all_outs.push(self.out_neighbors(into).as_ref().expect("into exists").clone());
        } else {
            return Err(ProcessingError::GraphError("Some node could not be removed.".to_owned()));
        }
        for f in from {
            if let Some((ins, outs)) = self.remove_node(f) {
                self.add_edges(ins.iter().map(|inc| (*inc, into)));
                self.add_edges(outs.iter().map(|outg| (into, *outg)));
                all_ins.push(ins);
                all_outs.push(outs);
            } else {
                return Err(ProcessingError::GraphError("Some node could not be removed.".to_owned()));
            }
        }
        self.remove_edge(&(into, into));
        Ok((all_ins, all_outs))
    }

    /// Contracts `node` adding edges between the incoming and the outgoing neighbors of `node`. 
    /// Returns the incoming neighbors, the outgoing neighbors and a set of already existing edges
    /// between outgoing and incoming neighbors, or `None` if `node` was deleted.
    ///
    /// # Panics 
    /// Panics if `node` is out of bounds.
    pub fn contract_node(&mut self, node: usize) -> Option<(FxHashSet<usize>, FxHashSet<usize>, HashSet<(usize, usize)>)> {
        if let Some(ins) = self.in_neighbors(node) {
            let mut doubles = HashSet::new();
            let ins_c = ins.clone();
            let outs_c = self.out_neighbors(node).as_ref().expect("`node` exists").clone();
            for trg in &outs_c {
                self.remove_edge(&(node, *trg));
            }
            for src in &ins_c {
                self.remove_edge(&(*src, node));
                for trg in &outs_c {
                    if !self.add_edge_checked((*src, *trg)) {
                        doubles.insert((*src,*trg));
                    }
                }
            }
            self.remove_node(node);
            return Some((ins_c, outs_c, doubles))
        }
        None
    }
    
}

// Reading graphs and building subgraphs
impl Digraph {
    pub fn read_graph<R: BufRead>(gr: R) -> Result<Self, ImportError> {
        let (lines, _comments): (Vec<_>, Vec<_>) = gr.lines()
            .partition(|l| {
                if let Ok(line) = l {
                    // separate empty lines and comment lines
                    !line.starts_with("% ")
                } else {
                    true
                }
            });
        let mut lines = lines.into_iter();
        let n = {
            let line = lines.next().ok_or(ImportError::InputMalformedError)??;
            let mut s = line.split(' ');
            let n: usize = s.next().ok_or(ImportError::InputMalformedError)?.parse()?;
            let _: usize = s.next().ok_or(ImportError::InputMalformedError)?.parse()?;
            if let Some("0") = s.next() {} else { return Err(ImportError::InputMalformedError); }
            if s.next().is_some() { return Err(ImportError::InputMalformedError); }
            n
        };
        let mut out_list = Vec::new();
        let mut in_list = vec![Some(FxHashSet::default()); n];
        let mut num_nodes = 0;
        for line in lines {
            let line = line?;
            let mut out_neighbors = FxHashSet::default();   
            if !line.is_empty() {
                let s = line.trim().split(' ');
                for unparsed_n in s {
                    let neighbor = unparsed_n.parse::<usize>()? - 1;
                    out_neighbors.insert(neighbor);
                    in_list[neighbor].as_mut().expect("No node ids over `n`.").insert(num_nodes);
                }
            }
            out_list.push(Some(out_neighbors));
            num_nodes += 1;
        }
        if num_nodes != n { return Err(ImportError::InputMalformedError); }

        Ok(Digraph {
            out_list,
            in_list,
        })
    }

    /// Write graph as metis file to a writer.
    pub fn write_graph<W: Write>(&self, mut out: W) -> Result<(), io::Error> {
        let mut i = 0;
        let node_name: Vec<_> = self.out_list.iter()
            .map(|n| {
                if n.is_some() {
                    i += 1;
                    i
                } else {
                    123321 // dummy value
                }
            })
            .collect();
        writeln!(out, "{:?} {:?} 0", self.num_nodes(), self.num_edges())?;
        for out_neighbors in self.out_list.iter().filter_map(|o| o.as_ref()) {
            let line = out_neighbors.iter().map(|nei| node_name[*nei].to_string()).collect::<Vec<_>>().join(" ");
            writeln!(out, "{}", line)?;
        }
        Ok(())
    }

    /// Write strong edge graph as .gr file to a writer.
    pub fn write_strong_graph_to_gr<W: Write>(&self, mut out: W) -> Result<(), io::Error> {
        writeln!(out, "p td {:?} {:?}", self.num_reserved_nodes(), self.strong_edges().count())?;
        for strong in self.strong_edges() {
            writeln!(out, "{:?} {:?}", strong.0, strong.1)?;
        }
        Ok(())
    }

    /// Transforms `other` (G(V,E)) to a graph G'(V',E') where each original node is replaced by an
    /// edge.
    /// Such that  V' = V U {v_n, v_(n+1),..., v_(n+n-1)}
    /// and E' = {(v_i, v_(i+n)) | i \in {0,...,n-1}} U {(v_(j+n), v_i) | (v_j, v_i) \in E} 
    pub (crate) fn transform_nodes_to_edges(other: &Self) -> Self {
        let num_old_nodes = other.num_nodes();
        let mut new_ins: Vec<_> = other.in_list.clone();
        new_ins.reserve_exact(num_old_nodes);
        let mut new_outs: Vec<Option<FxHashSet<usize>>> = Vec::with_capacity(num_old_nodes);
        let add_id = other.num_reserved_nodes();
        for i in other.nodes() {
            new_ins[i] = Some(other.in_list[i].as_ref().unwrap().iter().map(|node| node+add_id).collect::<FxHashSet<usize>>());
            new_ins.push(Some(vec![i].into_iter().collect::<FxHashSet<usize>>()));
            new_outs.push(Some(vec![i + num_old_nodes].into_iter().collect::<FxHashSet<usize>>()));
        }
        new_outs.append(other.out_list.clone().as_mut());
        Digraph {
            out_list: new_outs,
            in_list: new_ins,
        }
    }

    /// Transforms `other` (G(V,E)) to a graph G'(V',E') where the set of `to_transform` nodes are
    /// replaced edges.
    /// Such that  V' = V U {v_n, v_(n+1),..., v_(n+|`to_transform`|-1)}
    /// and E' = {(v_i, v_(i+n)) | i \in `to_transform`} U 
    ///          {(v_(j+n), v_i) | (v_j, v_i) \in E and j \in `to_transform`} 
    pub (crate) fn transform_some_nodes_to_edges(other: &Self, to_transform: &mut Vec<usize>) -> Self {
        to_transform.sort_unstable();
        let mut cur_index = other.num_reserved_nodes();
        let old_index = cur_index;
        let mut new_ins: Vec<_> = other.in_list.clone();
        new_ins.reserve_exact(to_transform.len());
        let mut new_outs: Vec<_> = other.out_list.clone();
        new_outs.reserve_exact(to_transform.len());
        let mut change_pairs: Vec<(usize, usize)> = Vec::new();
        let mut add_index: FxHashMap<usize, usize> = FxHashMap::default();
        for (id, node) in to_transform.iter().enumerate() {
            add_index.insert(*node, id);
            for out in other.out_list[*node].as_ref().unwrap() {
                change_pairs.push((*out, *node));
            }
            new_ins.push(Some(vec![*node].into_iter().collect::<FxHashSet<usize>>()));
            new_outs.push(other.out_list[*node].clone());
            new_outs[*node] = Some(vec![cur_index].into_iter().collect::<FxHashSet<usize>>());
            cur_index += 1;
        }
        change_pairs.sort_by_key(|(i,_)| *i);
        // There is probably a nicer way to do this:
        let mut current = if !change_pairs.is_empty() {change_pairs[0].0} else {0};
        let mut to_change: FxHashSet<usize> = FxHashSet::default();
        for (out, inn) in change_pairs {
            if out != current {
                new_ins[current] = Some(other.in_list[current].as_ref().unwrap()
                                    .iter()
                                    .map(|node| {
                                        if to_change.contains(node) {
                                            old_index+add_index.get(node).unwrap()
                                        } else {
                                            *node
                                        }
                                    }).collect::<FxHashSet<usize>>());
                current = out;
                to_change = FxHashSet::default();
            }
            to_change.insert(inn); 
        }
        if !to_change.is_empty() {
            new_ins[current] = Some(other.in_list[current].as_ref().unwrap()
                                .iter()
                                .map(|node| {
                                    if to_change.contains(node) {
                                        old_index+add_index.get(node).unwrap()
                                    } else {
                                        *node
                                    }
                                }).collect::<FxHashSet<usize>>());
        }
        Digraph {
            out_list: new_outs,
            in_list: new_ins,
        }
    }

    /// Build a subgraph given the set `node_set` of nodes the subgraph will contain.
    pub fn build_subgraph(&self, node_set: &FxHashSet<usize>) -> Self {
        let in_list = self.in_list.iter().cloned()
            .enumerate()
            .map(|(i, opt_list)| {
                if node_set.contains(&i) {
                    opt_list.map(|list| list.intersection(node_set).copied().collect::<FxHashSet<usize>>())
                }
                else {
                    None
                }
            }).collect::<Vec<Option<FxHashSet<usize>>>>();
        let out_list = self.out_list.iter().cloned()
            .enumerate()
            .map(|(i, opt_list)| {
                if node_set.contains(&i) {
                    opt_list.map(|list| list.intersection(node_set).copied().collect::<FxHashSet<usize>>())
                }
                else {
                    None
                }
            }).collect::<Vec<Option<FxHashSet<usize>>>>();
        Digraph{
            in_list,
            out_list,
        }
    }

    /// Removes all strong edges from the graph.
    pub fn remove_all_strong_edges(&mut self) {
        let nodes = self.nodes().collect::<Vec<_>>();
        for node in nodes {
            let ins = self.in_list[node].as_ref().expect("`node` exists");
            let outs = self.out_list[node].as_ref().expect("`node` exists");
            let intersect = ins.intersection(outs).cloned().collect::<Vec<_>>();
            for strong_neighbor in intersect {
                self.remove_edge(&(node,strong_neighbor));
                self.remove_edge(&(strong_neighbor,node));
            }
        }
    }

    /// Splits `self` into severel subgraphs, each containing exectly one of the connected
    /// components of `self`
    ///
    /// Should be more efficient than `self.split_into_connected_components()`
    pub fn split_into_connected_components_alt(&self) -> Vec<Self> {
        let big_n = self.num_reserved_nodes();
        let mut marked: FxHashSet<usize> = FxHashSet::default();
        let mut sub_graphs: Vec<Self> = Vec::new();
        for node in self.nodes() {
            if marked.contains(&node) {
                continue
            }
            let mut in_list: Vec<Option<FxHashSet<usize>>> = vec![None;big_n];
            let mut out_list: Vec<Option<FxHashSet<usize>>> = vec![None;big_n];
            for reachable_node in self.undirected_reachable(node) {
                marked.insert(reachable_node);
                in_list[reachable_node] = self.in_list[reachable_node].clone();
                out_list[reachable_node] = self.out_list[reachable_node].clone();
            }
            sub_graphs.push(Digraph {
                in_list,
                out_list,
            });
        }
        sub_graphs
    }

    /// Splits `self` into severel subgraphs, each containing exectly one of the connected
    /// components of `self`
    pub fn split_into_connected_components(&self) -> Vec<Self> {
        let mut marked: FxHashSet<usize> = FxHashSet::default();
        let mut sub_graphs: Vec<Self> = Vec::new();
        for node in self.nodes() {
            if marked.contains(&node) {
                continue
            }
            let reachable = self.undirected_reachable(node);
            reachable.iter().for_each(|node| {marked.insert(*node);});
            let mut new_sub_graph = self.clone();
            new_sub_graph
                .nodes()
                .collect::<FxHashSet<usize>>()
                .difference(&reachable)
                .for_each(|node| {
                    new_sub_graph.in_list[*node] = None;
                    new_sub_graph.out_list[*node] = None;
                });
            sub_graphs.push(new_sub_graph);
        }
        sub_graphs
    }

    /// DFS that puts all finished nodes (except the first) onto a stack.
    fn stack_dfs(&self, source: usize, marked: &mut NodeSet, stack: &mut Vec<usize>) {
        for neigh in self.out_neighbors(source).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)") {
            if marked.contains(neigh) {
                continue;
            }
            marked.insert(*neigh);
            self.stack_dfs(*neigh, marked, stack);
            stack.push(*neigh);
        }
    }


    /// Reversed DFS that puts all found nodes into a set.
    fn reverse_set_dfs(&self, source: usize, marked: &mut NodeSet, set: &mut FxHashSet<usize>) {
        for neigh in self.in_neighbors(source).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)") {
            if marked.contains(neigh) {
                continue;
            }
            marked.insert(*neigh);
            self.reverse_set_dfs(*neigh, marked, set);
            set.insert(*neigh);
        }
    }

    /// Finds the strongly connected components of `self`.
    pub fn find_strongly_connected_components(&self) -> Vec<FxHashSet<usize>> {
        let mut marked: NodeSet = NodeSet::new();
        let mut stack: Vec<usize> = Vec::new(); // For step 2
        // Step 1: (Recursively) start a dfs from each node in the graph and save nodes on `stack`
        // when ever their dfs is finished.
        for node in self.nodes() {
            if !marked.contains(&node) {
                marked.insert(node);
                self.stack_dfs(node, &mut marked, &mut stack);
                stack.push(node);
            }
        }
        let mut marked: NodeSet = NodeSet::new();
        let mut sccs: Vec<FxHashSet<usize>> = Vec::new();
        // Step 2: Start dfs in reversed direction from the nodes in `stack` (last in first out)
        // and add all traversed (not marked) nodes to the scc of current node from `stack`.
        while !stack.is_empty() {
            let node = stack.pop().expect("`stack` is not empty");
            if !marked.contains(&node) {
                marked.insert(node);
                let mut scc: FxHashSet<usize> = vec![node].into_iter().collect();
                self.reverse_set_dfs(node, &mut marked, &mut scc);
                sccs.push(scc);
            }
        }
        sccs
    }

    /// Finds the strongly connected components of `self`.
    pub fn find_strongly_connected_components_iter(&self) -> Vec<FxHashSet<usize>> {
        let mut marked: NodeSet = NodeSet::new();
        let mut stack: Vec<Im> = Vec::new(); // For step 2
        // Step 1: (Recursively) start a dfs from each node in the graph and save nodes on `stack`
        // when ever their respective dfs is finished.
        let mut queue: Vec<Im> = Vec::new();
        let mut nodes: Vec<Im> = self.nodes().map(|n| Im::Itm(n)).collect();
        queue.append(&mut nodes);
        while !queue.is_empty() {
            match queue.pop().expect("`queue` is not empty"){
                Im::Itm(node) => {
                    if !marked.contains(&node) {
                        queue.push(Im::Marker(node));
                        marked.insert(node);
                        let mut neighs: Vec<Im> = self.out_neighbors(node).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)").iter().copied().map(|n| Im::Itm(n)).collect();
                        queue.append(&mut neighs);
                    }
                },
                Im::Marker(node) => stack.push(Im::Itm(node)),
            }
        }
        let mut marked: NodeSet = NodeSet::new();
        let mut sccs: Vec<FxHashSet<usize>> = Vec::new();
        // Step 2: Start dfs in reversed direction from the nodes in `stack` (last in first out)
        // and add all traversed (not marked) nodes to the scc of current node from `stack`.
        let mut scc: FxHashSet<usize> = FxHashSet::default();
        while !stack.is_empty() {
            match stack.pop().expect("`queue` is not empty"){
                Im::Itm(node) => {
                    if !marked.contains(&node) {
                        marked.insert(node);
                        if scc.is_empty() {
                            stack.push(Im::Marker(node));
                        }
                        scc.insert(node);
                        let mut neighs: Vec<Im> = self.in_neighbors(node).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)").iter().copied().map(|n| Im::Itm(n)).collect();
                        stack.append(&mut neighs);
                    }
                },
                Im::Marker(_) => {
                    sccs.push(scc);
                    scc = FxHashSet::default();
                },
            }
        }
        sccs
    }

    /// DFS that ignores strong edges and that puts all finished nodes (except the first) onto a stack.
    fn weak_stack_dfs(&self, source: usize, marked: &mut NodeSet, stack: &mut Vec<usize>) {
        for neigh in self.weak_out_neighbors(source).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)") {
            if marked.contains(neigh) {
                continue;
            }
            marked.insert(*neigh);
            self.weak_stack_dfs(*neigh, marked, stack);
            stack.push(*neigh);
        }
    }

    /// Reversed DFS that ignores strong edges and that puts all found nodes into a set.
    fn weak_reverse_set_dfs(&self, source: usize, marked: &mut NodeSet, set: &mut FxHashSet<usize>) {
        for neigh in self.weak_in_neighbors(source).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)") {
            if marked.contains(neigh) {
                continue;
            }
            marked.insert(*neigh);
            self.weak_reverse_set_dfs(*neigh, marked, set);
            set.insert(*neigh);
        }
    }

    /// Finds the strongly connected components of `self` ignoring all strong edges.
    pub fn find_weak_strongly_connected_components(&self) -> Vec<FxHashSet<usize>> {
        let mut marked: NodeSet = NodeSet::new();
        let mut stack: Vec<usize> = Vec::new(); // For step 2
        // Step 1: (Recursively) start a dfs from each node in the graph and save nodes on `stack`
        // when ever their dfs is finished.
        for node in self.nodes() {
            if !marked.contains(&node) {
                marked.insert(node);
                self.weak_stack_dfs(node, &mut marked, &mut stack);
                stack.push(node);
            }
        }
        let mut marked: NodeSet = NodeSet::new();
        let mut sccs: Vec<FxHashSet<usize>> = Vec::new();
        // Step 2: Start dfs in reversed direction from the nodes in `stack` (last in first out)
        // and add all traversed (not marked) nodes to the scc of current node from `stack`.
        while !stack.is_empty() {
            let node = stack.pop().expect("`stack` is not empty");
            if !marked.contains(&node) {
                marked.insert(node);
                let mut scc: FxHashSet<usize> = vec![node].into_iter().collect();
                self.weak_reverse_set_dfs(node, &mut marked, &mut scc);
                sccs.push(scc);
            }
        }
        sccs
    }

    /// Finds the strongly connected components of `self`.
    pub fn find_weak_strongly_connected_components_iter(&self) -> Vec<FxHashSet<usize>> {
        let mut marked: NodeSet = NodeSet::new();
        let mut stack: Vec<Im> = Vec::new(); // For step 2
        // Step 1: (Recursively) start a dfs from each node in the graph and save nodes on `stack`
        // when ever their dfs is finished.
        let mut queue: Vec<Im> = Vec::new();
        let mut nodes: Vec<Im> = self.nodes().map(|n| Im::Itm(n)).collect();
        queue.append(&mut nodes);
        while !queue.is_empty() {
            match queue.pop().expect("`queue` is not empty"){
                Im::Itm(node) => {
                    if !marked.contains(&node) {
                        queue.push(Im::Marker(node));
                        marked.insert(node);
                        let mut neighs: Vec<Im> = self.weak_out_neighbors(node).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)").iter().copied().map(|n| Im::Itm(n)).collect();
                        queue.append(&mut neighs);
                    }
                },
                Im::Marker(node) => stack.push(Im::Itm(node)),
            }
        }
        let mut marked: NodeSet = NodeSet::new();
        let mut sccs: Vec<FxHashSet<usize>> = Vec::new();
        // Step 2: Start dfs in reversed direction from the nodes in `stack` (last in first out)
        // and add all traversed (not marked) nodes to the scc of current node from `stack`.
        let mut scc: FxHashSet<usize> = FxHashSet::default();
        while !stack.is_empty() {
            match stack.pop().expect("`queue` is not empty"){
                Im::Itm(node) => {
                    if !marked.contains(&node) {
                        marked.insert(node);
                        if scc.is_empty() {
                            stack.push(Im::Marker(node));
                        }
                        scc.insert(node);
                        let mut neighs: Vec<Im> = self.weak_in_neighbors(node).as_ref().expect("current is either a node in the graph or a neighbor of an existing node (which makes it a node of the graph)").iter().copied().map(|n| Im::Itm(n)).collect();
                        stack.append(&mut neighs);
                    }
                },
                Im::Marker(_) => {
                    sccs.push(scc);
                    scc = FxHashSet::default();
                },
            }
        }
        sccs
    }

    /// Greedily looks for a maximal clique: 
    /// Starting with the node with the highest strong degree and then adding node with the biggest
    /// neighborhood intersection.
    pub fn greedy_max_clique(&self) -> FxHashSet<usize> {
        let mut greedy_max_clique = FxHashSet::default();
        if let Some((node, mut strong_neighbors)) = self.get_max_strong_degree_node_and_neighbors() {
            greedy_max_clique.insert(node);
            while !strong_neighbors.is_empty() {
                let mut max_intersect = FxHashSet::default();
                let mut max_neighbor = None;
                for neigh in &strong_neighbors {
                    if let Some(shared_strong_neighbors) = self.strong_neighbors_in(*neigh, &strong_neighbors) {
                        if shared_strong_neighbors.len() > max_intersect.len() || max_neighbor.is_none() {
                            max_intersect = shared_strong_neighbors;
                            max_neighbor = Some(*neigh);
                        }
                    }
                }
                if let Some(neigh) = max_neighbor {
                    greedy_max_clique.insert(neigh);
                    strong_neighbors = max_intersect;
                } else {
                    break
                }
            }
        }
        greedy_max_clique
    }

    /// Returns an inclusion wise maximal double path of even length where the inner nodes are only
    /// connected to the predecessor and successor, or `None` if no such path exists.
    ///
    /// Double circles are only returned if they are of even length.
    pub fn inclusion_maximal_even_double_path(&self) -> Option<Vec<usize>> {
        let mut marked = FxHashSet::default();
        'outer: for node in self.nodes() {
            if marked.contains(&node) {
                continue
            }
            marked.insert(node);
            if let Some(neighs) = self.is_link_node(node) {
                let mut pre = node; 
                let mut right = neighs[0];
                let mut to_the_right = Vec::new();
                'right: loop {
                    to_the_right.push(right);
                    // if we found double circle
                    if right == node {
                        if to_the_right.len() % 2 == 0 {
                            return Some(to_the_right);
                        } else {
                            continue 'outer
                        }
                    }
                    marked.insert(right);
                    if let Some(neighs) = self.is_link_node(right) {
                        if pre == neighs[0] {
                            pre = right;
                            right = neighs[1];
                        } else {
                            pre = right;
                            right = neighs[0];
                        }
                    } else {
                        break 'right
                    }
                }
                let mut pre = node; 
                let mut left = neighs[1];
                let mut to_the_left = Vec::new();
                'left: loop {
                    to_the_left.push(left);
                    marked.insert(left);
                    if let Some(neighs) = self.is_link_node(left) {
                        if pre == neighs[0] {
                            pre = left;
                            left = neighs[1];
                        } else {
                            pre = left;
                            left = neighs[0];
                        }
                    } else {
                        break 'left
                    }
                }
                // check if path is even
                if (to_the_left.len() + to_the_right.len() + 1) % 2 == 0 {
                    to_the_left.reverse();
                    to_the_left.push(node);
                    to_the_left.append(&mut to_the_right);
                    return Some(to_the_left);
                }
            }
        }
        None
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn read_digraph_test() {
        let gr = Cursor::new("7 10 0\n3 4 5 6 7\n3 4 5 6 7\n\n\n\n\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
    }

    #[test]
    fn find_daisies_test() {
        let gr = Cursor::new("7 13 0\n2 3 4 7\n1 4 5\n2\n1 2\n2\n1\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert_eq!(g.get_max_strong_degree_node_and_neighbors(),
        Some((1,vec![0,3,4].into_iter().collect::<FxHashSet<usize>>())));
    }


    #[test]
    fn transform_digraph_test() {
        // TODO Test subgraph
        // TODO Test transformation and subgraphing on subgraph (or a graph with deleted nodes)
        let gr = Cursor::new("2 2 0\n2\n1\n");
        let gr_res = Cursor::new("4 4 0\n3\n4\n2\n1\n");
        let g = Digraph::read_graph(gr);
        let g_res = Digraph::read_graph(gr_res).unwrap();
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let g_new = Digraph::transform_nodes_to_edges(&mut g);
        assert_eq!(g_new, g_res);
        let gr = Cursor::new("3 4 0\n2\n1 3\n1\n");
        let gr_res = Cursor::new("5 6 0\n4\n5\n1\n2\n1 3\n");
        let g = Digraph::read_graph(gr);
        let g_res = Digraph::read_graph(gr_res).unwrap();
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let g_new = Digraph::transform_some_nodes_to_edges(&mut g, &mut vec![0usize,1]);
        assert_eq!(g_new, g_res);
    }

    #[test]
    fn find_components_test() {
        let gr = Cursor::new("10 10 0\n4\n9\n\n3 5\n1\n5 7\n8\n\n10\n2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let connected_comps = g.split_into_connected_components();
        assert_eq!(connected_comps.len(), 2);
        g.remove_node(1);
        g.remove_node(8);
        g.remove_node(9);
        let first_cp = &connected_comps[0];
        assert_eq!(first_cp, &g);
    }

    #[test]
    fn find_components_alt_test() {
        let gr = Cursor::new("10 10 0\n4\n9\n\n3 5\n1\n5 7\n8\n\n10\n2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let connected_comps = g.split_into_connected_components_alt();
        assert_eq!(connected_comps.len(), 2);
        g.remove_node(1);
        g.remove_node(8);
        g.remove_node(9);
        let first_cp = &connected_comps[0];
        assert_eq!(first_cp, &g);
    }

    #[test]
    fn greedy_max_clique_test() {
        let gr = Cursor::new("7 22 0\n2\n1 3 4 5\n2 4 5\n2 3 5 6 7\n2 3 4 7\n4 7\n4 6 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert_eq!(g.strong_neighbors_in(3, &vec![0,2,3,4].into_iter().collect::<FxHashSet<usize>>()).unwrap(), vec![2,4].into_iter().collect::<FxHashSet<usize>>());
        assert_eq!(g.greedy_max_clique(), vec![1,2,3,4].into_iter().collect::<FxHashSet<usize>>());
    }

    #[test]
    fn find_scc_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n\
                             7 8 9\n8 9\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let scc = g.find_strongly_connected_components();
        assert_eq!(scc.len(), 5);
        assert_eq!(scc[1],vec![1usize,2,3,4].into_iter().collect());
        assert_eq!(scc[4],vec![7usize,8,9,10].into_iter().collect());
    }

    #[test]
    fn find_scc_iter_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n\
                             7 8 9\n8 9\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let scc = g.find_strongly_connected_components_iter();
        assert_eq!(scc.len(), 5);
        assert_eq!(scc[1],vec![1usize,2,3,4].into_iter().collect());
        assert_eq!(scc[4],vec![7usize,8,9,10].into_iter().collect());
    }

    #[test]
    fn scc_graph_is_disconnected_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n\
                             7\n1\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert!(g.scc_graph_is_disconnected());
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n\
                             7 8\n1\n9\n11\n7 8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert!(!g.scc_graph_is_disconnected());
    }

    #[test]
    fn inclusion_maximal_even_double_path_test() {
        let gr = Cursor::new("9 17 0\n2\n1 3 9\n2 4\n3 5\n6 7\n4 5\n\
                            5 8\n7 9\n8\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let path = g.inclusion_maximal_even_double_path();
        assert!(path.is_some());
        let path = path.unwrap();
        assert_eq!(path.len(), 4);
        let gr = Cursor::new("9 18 0\n2 7\n1 3 9\n2 4\n3 5\n6 7\n4 5\n\
                            5 8\n7 9\n8\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let path = g.inclusion_maximal_even_double_path();
        assert!(path.is_none());
        let gr = Cursor::new("9 18 0\n2 5\n1 3\n2 4\n3 5\n4 1\n7 9\n\
                            8 6\n7 9\n9 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let path = g.inclusion_maximal_even_double_path();
        assert!(path.is_some());
        let path = path.unwrap();
        assert_eq!(path.len(), 4);
    }

    #[test]
    fn count_petals_test() {
        let gr = Cursor::new("20 26 0\n2 3 4 9 14\n1\n1\n5\n6 10\n7 15\n8\n\
                            1\n7\n11\n12\n13\n1\n6\n16\n17\n18\n19\n20\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert_eq!(g.count_petals(0), 5);
        let gr = Cursor::new("12 16 0\n2 5 12\n3\n4 6\n1\n8\n7\n1\n4 9\n\
                            10\n1\n12\n4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert_eq!(g.count_petals(0), 3);
    }

    #[test]
    fn count_petals_left_over_test() {
        let gr = Cursor::new("18 25 0\n2 3 5 6 8 16\n3\n1 4\n1\n7\n5\n6\n\
                            9\n10 12\n11\n1\n13\n14\n15\n1\n17\n18\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let (count, left_over) = g.count_petals_left_over(0);
        assert_eq!(count, 3);
        assert_eq!(left_over.num_nodes(), 5);
        assert_eq!(left_over.num_edges(), 3);
    }

    #[test]
    fn is_weak_edge_test() {
        let gr = Cursor::new("4 4 0\n2\n3\n2\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        assert!(g.is_weak_edge((0,1)));
        assert!(!g.is_weak_edge((1,0)));
        assert!(!g.is_weak_edge((1,2)));
        assert!(!g.is_weak_edge((2,3)));
    }

    #[test]
    fn find_split_set_test() {
        let gr = Cursor::new("14 22 0\n3 8\n1\n2 4\n2 6\n4 6\n\
                             7\n5\n9\n10 11\n8\n10 13\n11 13\n\
                             14\n7 12\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let split_set = g.find_split_set();
        assert_eq!(split_set.len(),2);
    }

}
