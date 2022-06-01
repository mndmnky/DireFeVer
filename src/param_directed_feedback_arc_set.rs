use crate::digraph::{Digraph};
use crate::param_directed_feedback_vertex_set::ParamDirectedFeedbackVertexSetInstance;
use crate::skew_edge_multicut::SkewEdgeMulticutInstance;
use std::convert::From;
use std::collections::HashSet;
use fxhash::{FxHashSet};
use itertools::Itertools;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ParamDirectedFeedbackArcSetInstance {
    pub (crate) graph: Digraph,
    pub (crate) node_edges: Option<usize>,
    pub (crate) k: usize,
}

impl ParamDirectedFeedbackArcSetInstance {

    /// Returns a new instance of the parametrized directed feedback arc set problem.
    pub fn new(graph: Digraph, k: usize) -> Self {
        ParamDirectedFeedbackArcSetInstance {
            graph,
            node_edges: None,
            k,
        }
    }

    /// Solves the directed feedback arc set problem by iterative compression. 
    ///
    /// Unclear whether this works probably on an instance that is not a converted directed
    /// feedback vertex set instance.
    pub fn solve(&self) -> Option<HashSet<(usize, usize)>> {
        // Find a proper starting subgraph.
        let mut compression_instance = self.build_init_subgraph();
        // If the subgraph is already complete and the heuristic solution is equal or lower than
        // `k`, we can already return the solution.
        if compression_instance.sub_result.len() <= self.k {
            if self.node_edges.is_some() {
                return Some(self.convert_vertex_to_arc_solution(&compression_instance.sub_result))
            } else {
                return compression_instance.solve()
            }
        }
        loop {
            // Solve the compression instance 
            if let Some(inter_solution) = compression_instance.solve() {
                // Convert intermediate arc solution to vertex solution.
                let vertex_sol: FxHashSet<usize> = inter_solution.iter().copied().map(|(_, trg)| trg).collect();
                if !self.extend_subgraph(&vertex_sol, &mut compression_instance) {
                    if self.node_edges.is_some(){
                        return Some(self.convert_vertex_to_arc_solution(&compression_instance.sub_result))
                    } else {
                        return compression_instance.solve();
                    }
                }
            } else {
                return None
            }
        }
    }

    /// Converts vertex solution of `self` into arc solution of `self`. This only works if `self`
    /// is a converted directed feedback vertex instance.
    ///
    /// # Panics
    /// Panics if `self` is no converted directed vertex feedback instance.
    fn convert_vertex_to_arc_solution(&self, vertex_sol: &FxHashSet<usize>) -> HashSet<(usize, usize)> {
        let threshold = self.node_edges.unwrap();
        vertex_sol.iter().map(|node| {
            if *node >= threshold {
                (node - threshold, *node)
            } else {
                (*node, node + threshold)
            }
        }).collect::<HashSet<(usize, usize)>>()
    }

    /// Extends the subgraph of `compression_instance` until the sub result is equal to `k + 1` or
    /// stops when no more nodes can be added. 
    ///
    /// Returns `false` if the sub result is smaller than `k + 1` but no more node can be added.
    fn extend_subgraph(&self,
                       intermediate_solution: &FxHashSet<usize>,
                       compression_instance: &mut ParamDirectedFeedbackArcSetCompressionInstance
                       ) -> bool { 
        // Get node_set and intermediate solution
        let old_node_set: FxHashSet<usize> = compression_instance.subgraph.nodes().collect();
        let mut new_node_set = old_node_set.clone();
        let mut sub_result = intermediate_solution.clone();
        // get unmarked_subnodes
        let mut unmarked_subnodes: FxHashSet<usize> = new_node_set.difference(intermediate_solution).copied().collect();
        let mut nodes = self.graph.nodes().filter(|n| !old_node_set.contains(n));
        let mut next_node = nodes.next();
        // add nodes as long as we have less or equal to `k` + 1 elements in W
        while next_node.is_some() && sub_result.len() <= self.k {
            let node = next_node.expect("Because we checked");
            new_node_set.insert(node);
            // If a node has a in neighbor and an out neighbor in `unmarked_subnodes` add it to
            // `sub_result` otherwise add it to `unmarked_subnodes`.
            if !self.graph.in_neighbors_in(node, &unmarked_subnodes).unwrap_or_default().is_empty() && !self.graph.out_neighbors_in(node, &unmarked_subnodes).unwrap_or_default().is_empty() {
                sub_result.insert(node);
            } else {
                unmarked_subnodes.insert(node);
            }
            next_node = nodes.next();
        }
        // build subgraph or extend if there is an efficient way to do so
        compression_instance.subgraph = self.graph.build_subgraph(&new_node_set);
        compression_instance.sub_result = sub_result;
        if compression_instance.sub_result.len() < self.k + 1 {
            return false
        }
        true
    }

    /// Heuristic to build initial subgraph for the compression instance.
    fn build_init_subgraph(&self) -> ParamDirectedFeedbackArcSetCompressionInstance {
        let mut subgraph_nodes: FxHashSet<usize> = FxHashSet::default();
        let mut sub_result: FxHashSet<usize> = FxHashSet::default();
        let mut nodes = self.graph.nodes();
        let mut next_node = nodes.next();
        let mut unmarked_subnodes: FxHashSet<usize> = FxHashSet::default();
        // Add nodes until there are no more nodes left or the size of `sub_result` is equal to
        // `k` + 1
        while next_node.is_some() && sub_result.len() <= self.k {
            let node = next_node.unwrap();
            subgraph_nodes.insert(node);
            // If a node has a in neighbor and an out neighbor in `unmarked_subnodes` add it to
            // `sub_result` otherwise add it to `unmarked_subnodes`.
            if !self.graph.in_neighbors_in(node, &unmarked_subnodes).unwrap_or_default().is_empty() && !self.graph.out_neighbors_in(node, &unmarked_subnodes).unwrap_or_default().is_empty() {
                sub_result.insert(node);
            } else {
                unmarked_subnodes.insert(node);
            }
            next_node = nodes.next();
        }
        let subgraph: Digraph = Digraph::build_subgraph(&self.graph, &subgraph_nodes);
        ParamDirectedFeedbackArcSetCompressionInstance::new(self.k, subgraph, sub_result)
    }

}

impl From<ParamDirectedFeedbackVertexSetInstance> for ParamDirectedFeedbackArcSetInstance {

    /// Converts an instance of the parametrized directed feedback vertex set problem to an
    /// instance of the parametrized directed feedback arc set problem.
    fn from(item: ParamDirectedFeedbackVertexSetInstance) -> Self {
        let graph = Digraph::transform_nodes_to_edges(&item.graph);
        let node_edges = Some(item.graph.num_nodes());
        let k = item.k;
        ParamDirectedFeedbackArcSetInstance{
            graph,
            node_edges,
            k,
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ParamDirectedFeedbackArcSetCompressionInstance {
    pub (crate) subgraph: Digraph,
    pub (crate) sub_result: FxHashSet<usize>,
    pub (crate) k: usize,
}

impl ParamDirectedFeedbackArcSetCompressionInstance {

    /// Returns an instance of the directed feedback arc set compression problem.
    pub fn new(k: usize, subgraph: Digraph, sub_result: FxHashSet<usize>) -> Self {
        ParamDirectedFeedbackArcSetCompressionInstance{
            subgraph,
            sub_result,
            k,
        }
    }

    /// Solves a directed feedback arc set compression instance by transforming it into a skew edge
    /// multicut instance.
    ///
    /// The result might be suboptimal if `self.sub_result.len()` < k + 1.
    /// Undefined behavior if `self.sub_result.len()` > k + 1.
    fn solve(&self) -> Option<HashSet<(usize, usize)>> {
        let skew_graph = Digraph::transform_some_nodes_to_edges(&self.subgraph, &mut self.sub_result.clone().into_iter().collect());
        let pairs: Vec<(usize, usize)> = self.sub_result.iter().map(|node| {
            (*skew_graph.out_neighbors(*node).as_ref().unwrap().iter().next().unwrap(),*node)
        }).collect();
        // permute `sub_result`:
        for perm in pairs.iter().permutations(self.k + 1) {
            // Transform subgraph to a skew edge multicut instance
            let node_pairs: Vec<(usize, usize)> = perm.into_iter().copied().collect();
            let mut skew_problem = SkewEdgeMulticutInstance::new(skew_graph.clone(), node_pairs, self.k);
            // Solve skew_edge_multicut instance:
            if let Some(solution) = skew_problem.solve() {
                return Some(solution)
            }
        }
        None
    }

}


#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use crate::digraph::Digraph;
    use super::*;

    #[test]
    fn solve_compression_instance_test() {
        let gr = Cursor::new("18 24 0\n10\n11\n12\n13\n14\n15\n16\n\
                                17\n18\n3\n1 4\n2 5\n3\n4 6\n\
                                5 7 9\n9\n6\n8 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let sk = Cursor::new("21 27 0\n10\n11\n12\n13\n14\n15\n16\n\
                                17\n18\n3\n19\n2 5\n3\n\
                                20\n5 7 9\n9\n6\n21\n1 4\n\
                                4 6\n8 6\n");
        let sk_c = Digraph::read_graph(sk);
        assert!(sk_c.is_ok());
        let sk_c = sk_c.unwrap();
        let skew_graph = Digraph::transform_some_nodes_to_edges(&g, &mut vec![10usize, 13, 17]);
        assert_eq!(sk_c, skew_graph);
        let dfasci = ParamDirectedFeedbackArcSetCompressionInstance::new(2, g, vec![10usize, 13, 17].into_iter().collect::<FxHashSet<usize>>());
        let solution = dfasci.solve();
        assert_eq!(solution, Some(vec![(2usize, 11usize), (5,14)].into_iter().collect::<HashSet<(usize, usize)>>()));
    }

    #[test]
    fn solve_param_directed_feedback_arc_instance_test() {
        let gr = Cursor::new("18 24 0\n10\n11\n12\n13\n14\n15\n16\n\
                                17\n18\n3\n1 4\n2 5\n3\n4 6\n\
                                5 7 9\n9\n6\n8 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let dfasi = ParamDirectedFeedbackArcSetInstance {
                graph: g,
                node_edges: Some(9),
                k: 2,
            };
        let solution = dfasi.solve();
        assert_eq!(solution, Some(vec![(2usize, 11usize), (5,14)].into_iter().collect::<HashSet<(usize, usize)>>()));
    }

}
