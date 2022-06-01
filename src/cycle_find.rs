use crate::digraph::Digraph;
use std::collections:: VecDeque;
use fxhash::FxHashSet;
use crate::other_ds::NodeSet;

impl Digraph{

    // TODO: Struggle with loops

    /// Finds a set of node disjoint cycles containing `node`. Always adding the smallest possible
    /// node disjoint cycle to the set.
    ///
    /// Returns a vector of cycles and the left over graph after deleting all those cycles.
    pub fn smallest_first_heuristic_node(&self, node: usize) -> Option<(Vec<Vec<usize>>, Self)> {
        if self.has_node(node) {
            let mut clone_graph = self.clone();
            let mut cycles = Vec::new();
            while clone_graph.reduce_to_cycles() && clone_graph.has_node(node) {
                if let Some(cycle) = clone_graph.smallest_simple_cycle(node) {
                    if cycle.is_empty() {
                        clone_graph.remove_edge(&(node, node)); // In case there are still loops left (which shouldn't be).
                    }
                    clone_graph.remove_nodes(cycle.clone());
                    cycles.push(cycle);
                } else {
                    break
                }
            }
            return Some((cycles, clone_graph))
        }
        None
    }

    //pub fn find_reducible_double_cycles_node(&self, node: usize) -> Option<HashSet<Vec<usize>>> {
    //    if !self.has_node(node) {
    //        return None
    //    }
    //    let mut used = FxHashSet::new();
    //    let mut visited = FxHashSet::new();
    //    // Let nodes with st_deg == 2 and wk_deg == 0 be called link nodes
    //}

    /// Finds the smallest simple cycle containing `node` if one exists. 
    pub fn smallest_simple_cycle(&self, node: usize) -> Option<Vec<usize>> {
        let mut marked: NodeSet = vec![node].into_iter().collect();
        let mut predecessors: Vec<Option<usize>> = vec![None; self.num_reserved_nodes()];
        let mut queue: VecDeque<usize> = vec![node].into_iter().collect();
        while !queue.is_empty() {
            let current = queue.pop_front().expect("`queue` is not empty");
            for neigh in self.out_neighbors(current).as_ref().expect("current is either `node` or a neighbor of an existing node") {
                if !marked.contains(neigh) {
                    queue.push_back(*neigh);
                    marked.insert(*neigh);
                    predecessors[*neigh] = Some(current);
                }
                if neigh == &node {
                    let mut cycle = Vec::new();
                    if current == node {
                        return Some(cycle) // In case of loop.
                    }
                    let mut tail = Some(current);
                    loop {
                        cycle.push(tail.expect("`tail` is some until `node` is found."));
                        tail = predecessors[tail.expect("`tail` is some until `node` is found.")];
                        if tail == Some(node) {
                            return Some(cycle)
                        }                         
                    }
                }
            }
        }
        None
    }


    /// Finds a set of node disjoint cycles containing `node`. Tries to find cycles wih nodes of
    /// small degree first.
    ///
    /// Returns a vector of cycles and the left over graph after deleting all those cycles.
    pub fn small_degree_heuristic_node(&self, node: usize) -> Option<(Vec<Vec<usize>>, Self)> {
        if self.has_node(node) {
            let mut clone_graph = self.clone();
            let mut cycles = Vec::new();
            while clone_graph.reduce_to_cycles() && clone_graph.has_node(node) {
                if let Some(cycle) = clone_graph.smallest_degree_simple_cycle_heuristic(node) {
                    if cycle.is_empty() {
                        clone_graph.remove_edge(&(node, node)); // In case there are still loops left (which shouldn't be).
                    }
                    clone_graph.remove_nodes(cycle.clone());
                    cycles.push(cycle);
                } else {
                    break
                }
            }
            return Some((cycles, clone_graph))
        }
        None
    }

    /// Tries to find the simple cycle containining nodes with low degree and `node` if one exists. 
    fn smallest_degree_simple_cycle_heuristic(&self, node: usize) -> Option<Vec<usize>> {
        let mut marked: NodeSet = vec![node].into_iter().collect();
        let mut predecessors: Vec<Option<usize>> = vec![None; self.num_reserved_nodes()];
        let mut queue: Vec<usize> = vec![node];
        while !queue.is_empty() {
            let current = queue.pop().expect("`queue` is not empty");
            // Get neighbors sorted by min max degree
            let mut outs: Vec<&usize> = self.out_neighbors(current).as_ref().expect("current is either `node` or a neighbor of an existing node").iter().collect();
            outs.sort_unstable_by_key(|out| self.degree(**out));
            outs.reverse();
            for neigh in outs{
                if !marked.contains(neigh) {
                    queue.push(*neigh);
                    marked.insert(*neigh);
                    predecessors[*neigh] = Some(current);
                }
                if neigh == &node {
                    let mut cycle = Vec::new();
                    if current == node {
                        return Some(cycle) // In case of loop.
                    }
                    let mut tail = Some(current);
                    loop {
                        cycle.push(tail.expect("`tail` is some until `node` is found."));
                        tail = predecessors[tail.expect("`tail` is some until `node` is found.")];
                        if tail == Some(node) {
                            return Some(cycle)
                        }                         
                    }
                }
            }
        }
        None
    }

    /// Reduces the graph to its cycles, by applying the simple rule that removes source and sinks,
    /// and splits the graph into strongly connected components.
    ///
    /// Returns true if there still exist cycles and false otherwise.
    ///
    /// TODO removal of sources and sinks should suffice
    pub fn reduce_to_cycles(&mut self) -> bool {
        let mut changed = true;
        while changed {
            changed = false;
            let nodes = self.nodes().collect::<Vec<usize>>();
            for node in nodes {
                // Remove source and sinks:
                if self.in_degree(node).expect("`node` is in graph.nodes()") == 0 || self.out_degree(node).expect("`node` is in graph.nodes()") == 0 {
                    self.remove_node(node);
                    changed = true;
                }
            }
        }
        if self.num_nodes() > 0 {
            // Split into strongly connected components:
            let sccs = self.find_strongly_connected_components();
            let matter_sccs: Vec<FxHashSet<usize>> = sccs.into_iter()
                .filter(|scc| {
                    if scc.len() == 1 {
                        let elem = scc.iter().next().expect("`scc` holds one element");
                        self.remove_node(*elem);
                        false
                    } else {
                        true
                    }
                }).collect();
            for i in 0..matter_sccs.len() {
                for j in (i+1)..matter_sccs.len(){
                    let edges_between: Vec<_> = self.edges_between(&matter_sccs[i], &matter_sccs[j]).collect();
                    if !edges_between.is_empty() {
                        self.remove_edges(edges_between);
                    }
                }
            }
            return true
        }
        false
    }

    /// Checks if `self` contains a cycle.
    pub fn has_cycle(&self) -> bool {
        let mut graph_clone = self.clone();
        'outer: loop {
            let nodes = graph_clone.nodes().collect::<Vec<usize>>();
            for node in nodes {
                // Remove source and sinks:
                if graph_clone.in_degree(node).expect("`node` is in graph.nodes()") == 0 || graph_clone.out_degree(node).expect("`node` is in graph.nodes()") == 0 {
                    graph_clone.remove_node(node);
                    continue 'outer 
                }
            }
            break
        }
        !(graph_clone.num_nodes() == 0)
    }

    /// Returns a vector of the smallest cycle in reversed ordering, or a cycle of size at most 2 in `self`, or None if no cycle
    /// exists. 
    pub fn find_smallest_cycle(&self) -> Option<Vec<usize>> {
        let mut smallest = None;
        let mut size = 0;
        for node in self.nodes() {
            if let Some(mut cycle) = self.smallest_simple_cycle(node){
                if cycle.len() < 2 {
                    cycle.push(node);
                    return Some(cycle)
                }
                if smallest.is_none() || size > cycle.len() + 1 {
                    cycle.push(node);
                    size = cycle.len();
                    smallest = Some(cycle);
                }
            }
        }
        smallest
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn smallest_first_test() {
        let gr = Cursor::new("13 17 0\n2 5 7\n3\n4\n1 12\n6\n\
        7\n8\n9\n10\n11\n12\n1 13\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let (cycles, _left_overs) = g.smallest_first_heuristic_node(0usize).unwrap();
        assert_eq!(cycles, vec![vec![3usize, 2, 1], vec![11,10,9,8,7,6]]);
    }

    #[test]
    fn small_degree_first_test() {
        let gr = Cursor::new("6 9 0\n2 3\n3 4\n4 6\n5\n6\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let (cycles, _left_overs) = g.smallest_first_heuristic_node(0usize).unwrap();
        let (cycles2, _left_overs) = g.small_degree_heuristic_node(0usize).unwrap();
        assert_eq!(cycles, vec![vec![5usize, 2]]);
        assert_eq!(cycles2, vec![vec![5usize, 4, 3, 1]]);
    }

    #[test]
    fn smallest_or_3_cycle_test() {
        let gr = Cursor::new("7 9 0\n2\n3\n4 7\n5\n4 6\n1\n2\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let cycle = g.find_smallest_cycle();
        assert_eq!(cycle, Some(vec![4usize, 3]));
    }

}
