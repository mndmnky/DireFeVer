use crate::digraph::Digraph;
use std::collections::{VecDeque, HashSet, HashMap};
use fxhash::FxHashSet;
use crate::other_ds::NodeSet;
use crate::digraph::Im;

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
            clone_graph.reduce_to_cycles();
            while clone_graph.has_node(node) {
                if let Some(cycle) = clone_graph.smallest_simple_cycle(node) {
                    if cycle.is_empty() {
                        clone_graph.remove_edge(&(node, node)); // In case there are still loops left (which shouldn't be).
                    }
                    clone_graph.remove_nodes(cycle.clone());
                    cycles.push(cycle);
                } else {
                    break
                }
                clone_graph.reduce_to_cycles();
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
            clone_graph.reduce_to_cycles();
            while clone_graph.has_node(node) {
                if let Some(cycle) = clone_graph.smallest_degree_simple_cycle_heuristic(node) {
                    if cycle.is_empty() {
                        clone_graph.remove_edge(&(node, node)); // In case there are still loops left (which shouldn't be).
                    }
                    clone_graph.remove_nodes(cycle.clone());
                    cycles.push(cycle);
                } else {
                    break
                }
                clone_graph.reduce_to_cycles();
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

    /// Finds semi node disjoint cycles in `self` so that no node is incident to more than
    /// `max_count` cycles.
    pub fn find_semi_disjoint_cycles(&self, max_count: usize) -> Vec<Vec<usize>> {
        let mut cycles = Vec::new();
        let mut node_cycle_counter: Vec<Option<usize>> = (0..self.num_reserved_nodes())
            .map(|node| {
                if self.has_node(node) {
                    Some(0)
                } else {
                    None
                }
            }).collect();
        let mut graph_clone = self.clone();
        loop {
            if let Some((next,amount)) = node_cycle_counter.iter().enumerate().filter(|(_,n)| n.is_some()).map(|(i,n)|(i,n.expect("is some"))).min_by_key(|(_,n)| *n) {
                if amount >= max_count {
                    break;
                }
                if let Some(next_cycle) = graph_clone.find_semi_disjoint_cycle(next, max_count, &node_cycle_counter) {
                    // remove each edge in cycle from graph_clone
                    let start = next_cycle[0];
                    let mut pref = start;
                    if node_cycle_counter[start].filter(|c| c<&(max_count-1)).is_some() {
                        node_cycle_counter[start] = node_cycle_counter[start].map(|c| c+1);
                    } else {
                        node_cycle_counter[start] = None;
                    }
                    for node in &next_cycle[1..next_cycle.len()] {
                        assert!(graph_clone.remove_edge(&(pref,*node)));
                        if node_cycle_counter[*node].filter(|c| c<&(max_count-1)).is_some() {
                            node_cycle_counter[*node] = node_cycle_counter[*node].map(|c| c+1);
                        } else {
                            node_cycle_counter[*node] = None;
                        }
                        pref = *node;
                    }
                    assert!(graph_clone.remove_edge(&(pref,start)));
                    graph_clone.reduce_to_cycles_fix_list(&mut node_cycle_counter);
                    cycles.push(next_cycle);
                } else {
                    // If no cycle with `next` exists 
                    graph_clone.remove_node(next);
                    node_cycle_counter[next] = None;
                    graph_clone.reduce_to_cycles_fix_list(&mut node_cycle_counter);
                }
            } else {
                break;
            }
        }
        cycles
    }



    /// Finds a cycle in `self` containing `start_node` that uses only nodes that have been used less then `max_count`
    /// times.
    /// 
    /// Cycles are found with a DFS, visiting nodes that have been used less first.
    ///
    /// Returns a cycle containing `start_node` or `None` if no such cycle could have been found.
    pub fn find_semi_disjoint_cycle(&mut self, start_node: usize, max_count: usize, node_cycle_counter: &Vec<Option<usize>>) 
        -> Option<Vec<usize>> {
        let mut marked: NodeSet = NodeSet::new();
        let mut stack: Vec<usize> = Vec::new(); 
        let mut queue: Vec<Im> = Vec::new();
        let mut found_cycle = None;
        queue.push(Im::Itm(start_node));
        while !queue.is_empty() {
            match queue.pop().expect("`queue` is not empty"){
                Im::Itm(node) => {
                    if !marked.contains(&node) {
                        queue.push(Im::Marker(node));
                        marked.insert(node);
                        stack.push(node);
                        let outs_n = self.out_neighbors(node).as_ref().expect("`next` exists").clone();
                        let mut outs = outs_n.iter().filter_map(|neigh| {
                            if let Some(count) = node_cycle_counter[*neigh]{
                                if count < max_count {
                                    return Some((count, neigh))
                                }
                            }
                            None
                        }).collect::<Vec<_>>();
                        outs.sort_by_key(|(count,_)| *count);
                        //outs.reverse();
                        while !outs.is_empty() {
                            let neigh = outs.pop().expect("`outs` is not empty").1;
                            queue.push(Im::Itm(*neigh));
                        }
                    } else {
                        if node == start_node {
                            found_cycle = Some(stack);
                            break;
                        }
                    }
                },
                Im::Marker(node) => {
                    assert_eq!(stack.pop(), Some(node));
                },
            }
        }
        found_cycle
    }

    /// Greedily finds weak, node disjoint cycles of size `size` or smaller and returns them.
    pub fn find_disjoint_cycles_of_at_most_size(&self, size: usize) -> HashSet<Vec<usize>> {
        let mut graph_clone = self.clone();
        let mut cycles = HashSet::new();
        while graph_clone.num_nodes() != 0 {
            // get any node 
            let node = graph_clone.nodes().next().expect("`num_nodes()`>0");
            // do bfs into depth `size` 
            if let Some(cycle) = graph_clone.find_weak_cycle_of_size_or_lower(node, size) {
                cycles.insert(cycle.clone());
                graph_clone.remove_nodes(cycle.into_iter());
                graph_clone.reduce_to_cycles();
            } else {
                graph_clone.remove_node(node);
                graph_clone.reduce_to_cycles();
            }
        }
        return cycles;
    }

    /// Greedily finds a weak, node disjoint cycles of size `size` or smaller 
    /// containing `node` and returns it if one was found.
    pub fn find_weak_cycle_of_size_or_lower(&self, node: usize, size: usize) -> Option<Vec<usize>> {
        let mut dists: HashMap<usize,usize> = HashMap::new();
        let mut pred: HashMap<usize,usize> = HashMap::new();
        let mut queue: VecDeque<(usize, usize)> = VecDeque::new();
        let mut cycle = Vec::new();
        queue.push_front((node, node));
        pred.insert(node, node);
        //dists.insert(node, 0);
        while !queue.is_empty() {
            let (next, pref) = queue.pop_front().expect("is not empty");
            if next == node && pref != node {
                let mut pn = pref;
                cycle.push(pn);
                while pn != node {
                    pn = *pred.get(&pn).expect("`pn` was a predecessor");
                    cycle.push(pn);
                }
                cycle.reverse();
                return Some(cycle);
            }
            if let Some(_) = dists.get(&next) {
                continue;
            } else {
                if next == node {
                    dists.insert(next, 0);
                } else {
                    dists.insert(next, dists.get(&pref).expect("pushed `next`") + 1);
                }
                if dists.get(&next).expect("just set") == &size {
                    continue;
                }
                pred.insert(next, pref);
                queue.extend(self.weak_out_neighbors(next).expect("`next` exists").iter().map(|o| (*o, next)));
            }
        }
        return None
    }

    /// Greedily finds node disjoint transitive edge structures of the form (a,b) (a,c) (c,b).
    pub fn find_disjoint_transitive_structures(&self) -> Vec<(usize, usize, usize)> {
        let mut graph_clone = self.clone();
        let mut structs = Vec::new();
        while graph_clone.num_edges() > 2 {
            let mut found = None;
            let w_edges: Vec<_> = graph_clone.weak_edges().collect();
            for w_edge in w_edges {
                if let Some(tst) = graph_clone.find_transitive_structure_for_edge(w_edge) {
                    found = Some(tst);
                    break;
                } else {
                    graph_clone.remove_edge(&w_edge);
                }
            }
            if let Some(tst) = found {
                structs.push(tst);
                graph_clone.remove_node(tst.0);
                graph_clone.remove_node(tst.1);
                graph_clone.remove_node(tst.2);
            } else {
                return structs;
            }
        }
        structs
    }

    /// Find a transitive edge structure ((a,b) (a,c) (c,b)) with edges not in PIE such that `edge`
    /// is part of that structure.
    fn find_transitive_structure_for_edge(&self, edge: (usize, usize)) -> Option<(usize, usize, usize)> {
        assert!(self.edge_exists(&edge));
        assert!(!self.edge_exists(&(edge.1, edge.0)));
        let outs0 = self.weak_out_neighbors(edge.0).expect("`edge` exists");
        let outs1 = self.weak_out_neighbors(edge.1).expect("`edge` exists");
        let ins0 = self.weak_in_neighbors(edge.0).expect("`edge` exists");
        let ins1 = self.weak_in_neighbors(edge.1).expect("`edge` exists");
        let candidates = outs0.intersection(&ins1);
        if let Some(choice) = candidates.min_by_key(|c| self.degree(**c).expect("`c` exists")) {
            return Some((edge.0, edge.1, *choice))
        } else {
            let candidates = outs0.intersection(&outs1);
            if let Some(choice) = candidates.min_by_key(|c| self.degree(**c).expect("`c` exists")) {
                return Some((edge.0, *choice, edge.1))
            } else {
                let candidates = ins0.intersection(&ins1);
                if let Some(choice) = candidates.min_by_key(|c| self.degree(**c).expect("`c` exists")) {
                    return Some((*choice, edge.1, edge.0))
                }
            }
        }
        None
    }

    /// Reduces the graph to its cycles, by applying the simple rule that removes source and sinks.
    /// Writes `None` to `list_to_fix` if a node was removed.
    ///
    /// Returns true if there still exist cycles and false otherwise.
    pub fn reduce_to_cycles_fix_list(&mut self, list_to_fix: &mut Vec<Option<usize>>) -> bool {
        let mut changed = true;
        while changed {
            changed = false;
            let nodes = self.nodes().collect::<Vec<usize>>();
            for node in nodes {
                // Remove source and sinks:
                if self.in_degree(node).expect("`node` is in graph.nodes()") == 0 || self.out_degree(node).expect("`node` is in graph.nodes()") == 0 {
                    self.remove_node(node);
                    list_to_fix[node] = None;
                    changed = true;
                }
            }
        }
        false
    }

    /// Reduces the graph to its cycles, by applying the simple rule that removes source and sinks.
    ///
    /// Returns true if there still exist cycles and false otherwise.
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
        false
    }

    /// Reduces the graph to its cycles, by applying the simple rule that removes source and sinks,
    /// and splits the graph into strongly connected components.
    ///
    /// Returns true if there still exist cycles and false otherwise.
    pub fn reduce_to_sccs(&mut self) -> Option<Vec<FxHashSet<usize>>> {
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
            let sccs = self.find_strongly_connected_components_iter();
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
            return Some(matter_sccs)
        }
        None
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

    #[test]
    fn find_semi_disjoint_cycles_test() {
        let gr = Cursor::new("5 7 0\n2 4\n3 5\n1\n2\n1\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let cycles = g.find_semi_disjoint_cycles(2);
        dbg!(&cycles);
        assert_eq!(cycles.len(),2);
        let gr = Cursor::new("8 16 0\n2 5\n3 8\n4 7\n1 6\n4 8\n3 5\n2 6\n1 7\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let cycles = g.find_semi_disjoint_cycles(2);
        assert!((cycles.len() >= 4) && (cycles.len() <= 5));
    }

    #[test]
    fn find_weak_cycle_of_size_or_lower_test() {
        let gr = Cursor::new("7 9 0\n2\n3 4\n1\n5 6\n2\n7\n4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let cycles = g.find_disjoint_cycles_of_at_most_size(3);
        assert_eq!(cycles.len(),2);
    }

}
