use crate::digraph::{Digraph};
use fxhash::{FxHashSet};
use std::collections::{VecDeque, HashSet};

#[derive(Debug, Default, Clone, Eq, PartialEq)]
pub struct SkewEdgeMulticutInstance {
    graph: Digraph,
    ordered_node_pairs: Vec<(usize, usize)>,
    k: usize,
}

impl SkewEdgeMulticutInstance {

    /// Returns a instance of the skew edge multicut problem.
    pub fn new(graph: Digraph, ordered_node_pairs: Vec<(usize, usize)>, k: usize) -> Self {
        SkewEdgeMulticutInstance {
            graph,
            ordered_node_pairs,
            k,
        }
    }

    /// Reduces `self.k` by `amount`.
    fn reduce_k(&mut self, amount: usize) {
        self.k -= amount;
    }

    /// Increases `self.k` by `amount`.
    fn recover_k(&mut self, amount: usize) {
        self.k += amount;
    }

    /// Solve parametrized skew edge multicut instance.
    pub fn solve(&mut self) -> Option<HashSet<(usize, usize)>> {
        if self.ordered_node_pairs.is_empty() {
            return Some(HashSet::new())
        }
        let mut terminals: FxHashSet<usize> = self.ordered_node_pairs.iter().map(|(_,t)| t).copied().collect();
        let ordered_pais_copy = self.ordered_node_pairs.clone();
        let mut s = self.ordered_node_pairs.pop().unwrap();
        let mut s_set: FxHashSet<usize> = vec![s.0].into_iter().collect();
        while !self.graph.is_reachable_from(&s_set, &terminals){
            if self.ordered_node_pairs.is_empty() {
                return Some(HashSet::new())
            }
            terminals = self.ordered_node_pairs.iter().map(|(_,t)| t).copied().collect();
            s = self.ordered_node_pairs.pop().unwrap();
            s_set = vec![s.0].into_iter().collect();
        }
        // Get r_max 
        if let Some(mut r_max) = self.graph.find_r_set_max_with_k(&s_set, &terminals, self.k){
            // Enumerate important cuts
            let important_cuts = self.graph.enumerate_important_cuts(&mut r_max, &terminals, &mut HashSet::new(), self.k);
            // Branch on important cuts
            for cut in important_cuts {
                self.graph.remove_edges(cut.clone());
                self.reduce_k(cut.len());
                if let Some(mut sol) = self.solve() {
                    sol.extend(&cut);
                    return Some(sol)
                } else {
                    self.recover_k(cut.len());
                    self.graph.add_edges(cut);
                }
            }
        }
        self.ordered_node_pairs = ordered_pais_copy;
        None
    }
}

impl Digraph {

    /// Returns a vector with all important cuts.
    pub fn enumerate_important_cuts(&mut self, r_max: &mut FxHashSet<usize>, y_set: &FxHashSet<usize>, z_set: &mut HashSet<(usize,usize)>, k: usize) -> Vec<HashSet<(usize, usize)>> {
        let mut important_cuts: Vec<HashSet<(usize, usize)>> = Vec::new();
        let delta_x_set = self.get_cut(r_max);
        if delta_x_set.is_empty() {
            important_cuts.push(z_set.clone());
            return important_cuts
        }
        let edge = delta_x_set.into_iter().next().unwrap();
        if !y_set.contains(&edge.1){
            r_max.insert(edge.1);
            if let Some(mut new_r_max) = self.find_r_set_max_with_k(r_max, y_set, k) {
                // Check if no edge in `z_set` has both endpoints in `r_max`
                // Predict important cut
                if !z_set.iter().any(|z_edge| new_r_max.contains(&z_edge.1)) && 
                    self.predict_important_cut(&new_r_max, y_set, z_set, k) {
                    // branch
                    important_cuts.append(&mut self.enumerate_important_cuts(&mut new_r_max, y_set, z_set, k));
                }
            }
            r_max.remove(&edge.1);
        }
        z_set.insert(edge);
        // TODO: Do we need to recompute here?
        self.remove_edge(&edge);
        if let Some(mut new_r_max) = self.find_r_set_max_with_k(r_max, y_set, k) {
            // Check if no edge in Z has both endpoints in r_max
            // Predict important cut
            if !z_set.iter().any(|z_edge| new_r_max.contains(&z_edge.1)) && 
                self.predict_important_cut(&new_r_max, y_set, z_set, k) {
                // branch
                important_cuts.append(&mut self.enumerate_important_cuts(&mut new_r_max, y_set, z_set, k));
            }
        }
        self.add_edge(edge);
        z_set.remove(&edge);
        important_cuts
    }

    /// TODO: node set for multi graph (edge set?)
    /// Z is subset of DeltaX
    /// Lemma 8.15.
    pub fn predict_important_cut(&mut self, r_prime_max: &FxHashSet<usize>, y_set: &FxHashSet<usize>, z_set: &HashSet<(usize, usize)>, k: usize) -> bool {
        // TODO: could be avoided if we carry an original copy of the graph
        let mut g_orig = self.clone();
        g_orig.add_edges(z_set.clone());
        // Get delta r_prime_max
        let delta_r_prime_max: HashSet<(usize, usize)> = self.get_cut(r_prime_max);
        // Check if `r_prime_max` union `z_set` is an important cut of size at most `k` 
        // Check if `r_prime_max` union `z_set` is of size at most `k`
        if z_set.len() + delta_r_prime_max.len() > k {return false}
        // Check if `r_max` of a (`r_prime_max` union `z_set`) - `y_set` cut is equal to (`r_prime_max` union `z_set`)
        // r of `delta_r_prime_max` union `z_set` is simply `r_prime_max` since `z_set` is subset of delta `x_set`
        if let Some(r_check) = g_orig.find_r_set_max_with_k(r_prime_max, y_set, k) {
            // can be simply checked by their sizes
            return r_check.len() == r_prime_max.len()
        }
        false
    }

    /// Returns all edges (v, w) with v in `r_set` and w not in `r_set`.
    pub (crate) fn get_cut(&mut self, r_set: &FxHashSet<usize>) -> HashSet<(usize, usize)> {
        let mut cut: HashSet<(usize, usize)> = HashSet::new();
        for node in r_set {
            if let Some(neighborhood) = self.out_neighbors(*node) {
                cut.extend(neighborhood.iter().filter_map(|trg| {
                    if !r_set.contains(trg) {
                        Some((*node, *trg))
                    } else {
                        None
                    }
                }).collect::<HashSet<(usize,usize)>>());
            }
        }
        cut
    }

    /// Finds the node maximal set of all nodes that can not reach `y_set` after a minimum cut of size at
    /// most `k`.
    ///
    /// Attention! This does also add all other vertices that do not reach `y_set` after the cut, even if
    /// they are not reachable from `x_set`. (So far this makes no difference though)
    pub (crate) fn find_r_set_max_with_k(&self, x_set: &FxHashSet<usize>, y_set: &FxHashSet<usize>, k: usize) ->
        Option<FxHashSet<usize>> {
        // Start with `self` as residual network:
        let mut residual: Digraph = self.clone();
        let mut i = 0;
        // do BFS to find shortest augmenting path
        while let Some(path) = residual.find_augmenting_path(x_set, y_set) {
            if i == k {
                return None
            }
            i += 1; 
            // fix `residual`
            residual.rev_edges(path);
        }
        let rev_reachable: FxHashSet<usize> = residual.nodes_rev_reachable_from(y_set);
        let node_set: FxHashSet<usize> = self.nodes().collect();
        let r_max: FxHashSet<usize> = node_set.difference(&rev_reachable).copied().collect();
        Some(r_max)
    }

    /// Checks if `trg` is reachable from `src` in `self`.
    fn is_reachable_from(&self, src: &FxHashSet<usize>, trg: &FxHashSet<usize>) -> bool {
        let mut queue: VecDeque<usize> = VecDeque::new();
        let mut marked: FxHashSet<usize> = FxHashSet::default();
        for s in src {
            queue.push_back(*s);
            marked.insert(*s);
        }
        while !queue.is_empty(){
            let v: usize = queue.pop_front().unwrap();
            for node in self.out_neighbors(v).as_ref().unwrap().iter() {
                if !marked.contains(node) {
                    queue.push_back(*node);
                    marked.insert(*node);
                    if trg.contains(node) {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Returns a path from `src` to `trg` in `self`.
    pub (crate) fn find_augmenting_path(&mut self, src: &FxHashSet<usize>, trg: &FxHashSet<usize>) -> 
        Option<Vec<(usize, usize)>>{
        // TODO: Could be replaced by a `NodeSet` kind of structure.
        let mut predecessors: Vec<Option<usize>> = vec![None; self.num_reserved_nodes()];
        let mut queue: VecDeque<usize> = VecDeque::new();
        let mut path: Option<Vec<(usize, usize)>> = None;
        for s in src {
            queue.push_back(*s);
            predecessors[*s] = Some(*s);
        }
        while !queue.is_empty(){
            let v: usize = queue.pop_front().unwrap();
            self.out_neighbors(v).as_ref().unwrap().iter().for_each(|node| {
                if predecessors[*node].is_none() {
                    queue.push_back(*node);
                    predecessors[*node] = Some(v);
                    if trg.contains(node) {
                        let mut p: Vec<(usize, usize)> = Vec::new();
                        let mut t = *node;
                        loop {
                            if src.contains(&t) {break;}
                            let s = predecessors[t].unwrap();
                            p.push((s, t));
                            t = s;
                        }
                        path = Some(p);
                    }
                }
            });
            if path.is_some() {break;}
        }
        path
    }

    /// Returns all nodes that are reachable from `trg` by only using backward edges.
    pub (crate) fn nodes_rev_reachable_from(&mut self, trg: &FxHashSet<usize>) -> FxHashSet<usize> {
        let mut rev_reachable: FxHashSet<usize> = trg.clone();
        // Add `trg` to `queue`
        let mut queue: VecDeque<usize> = VecDeque::new();
        trg.iter().for_each(|node| {
            queue.push_back(*node);
        });
        // Until `queue` is empty:
        while !queue.is_empty(){
            // pop the front item of `queue`:
            let v = queue.pop_front().unwrap();
            // add `v`s incoming neighbors to `rev_reachable` and `queue` if
            // not marked:
            self.in_neighbors(v).as_ref().unwrap().iter().for_each(|node| {
                if !rev_reachable.contains(node) {
                    rev_reachable.insert(*node);
                    queue.push_back(*node);
                }
            });
        }
        rev_reachable
    }
}

#[cfg(test)]
mod tests {
    use crate::skew_edge_multicut::*;
    use std::io::Cursor;
    use crate::digraph::Digraph;
    use fxhash::{FxHashSet};

    #[test]
    fn min_cut_test() {
        let gr = Cursor::new("6 7 0\n2 3\n4 5\n4\n6\n6\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut x_set: FxHashSet<usize> = FxHashSet::default();
        let mut y_set: FxHashSet<usize> = FxHashSet::default();
        let res: FxHashSet<usize> = vec![0usize,1,2,3,4].into_iter().collect();
        x_set.insert(0);
        y_set.insert(5);
        assert!(g.find_r_set_max_with_k(&x_set, &y_set, 1).is_none());
        let opt_r_max = g.find_r_set_max_with_k(&x_set, &y_set, 2);
        assert_eq!(opt_r_max, Some(res));
        let gr = Cursor::new("9 9 0\n4 7\n4 5\n4\n6 9\n\n7 8\n\n\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let x: FxHashSet<usize> = vec![0usize,1,2].into_iter().collect();
        let y: FxHashSet<usize> = vec![6usize,7,8].into_iter().collect();
        let res: FxHashSet<usize> = vec![0usize,1,2,3,4].into_iter().collect();
        assert!(g.find_r_set_max_with_k(&x, &y, 2).is_none());
        let opt_r_max = g.find_r_set_max_with_k(&x, &y, 3);
        assert_eq!(opt_r_max, Some(res.clone()));
        let opt_r_max = g.find_r_set_max_with_k(&x, &y, 4);
        assert_eq!(opt_r_max, Some(res));
    }

    #[test]
    fn predict_important_cuts_test() {
        let gr = Cursor::new("19 26 0\n4\n5\n6\n7 8\n7 8\n8 9\n\
                             11\n10 12\n6\n13 14\n14\n14 15\n\
                             4 7 16\n16 17\n17 18 19\n\n\n\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let r_prime_max: FxHashSet<usize> = vec![0usize, 1, 2,3,4,5,6,7,8,10].into_iter().collect();
        let y: FxHashSet<usize> = vec![15usize, 16, 17, 18].into_iter().collect();
        let z: HashSet<(usize, usize)> = vec![(7usize, 9usize)].into_iter().collect();
        g.remove_edges(z.clone());
        assert!(g.predict_important_cut(&r_prime_max, &y, &z, 3));
    }

    #[test]
    fn important_cuts_test() {
        let gr = Cursor::new("19 26 0\n4\n5\n6\n7 8\n7 8\n8 9\n\
                             11\n10 12\n6\n13 14\n14\n14 15\n\
                             4 7 16\n16 17\n17 18 19\n\n\n\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let mut g = g.unwrap();
        let x: FxHashSet<usize> = vec![0usize, 1, 2].into_iter().collect();
        let y: FxHashSet<usize> = vec![15usize, 16, 17, 18].into_iter().collect();
        let r_max = g.find_r_set_max_with_k(&x, &y, 2);
        assert!(r_max.is_none());
        let mut r_max = g.find_r_set_max_with_k(&x, &y, 3).unwrap();
        let cuts = g.enumerate_important_cuts(&mut r_max, &y, &mut HashSet::new(), 3);
        let res: HashSet<(usize, usize)> = vec![(7usize, 9usize), (10,13), (7,11)].into_iter().collect();
        assert_eq!(cuts, vec![res.clone()]);
        let mut r_max = g.find_r_set_max_with_k(&x, &y, 6).unwrap();
        let cuts = g.enumerate_important_cuts(&mut r_max, &y, &mut HashSet::new(), 6);
        let res2: HashSet<(usize, usize)> = vec![(13usize, 15usize), (13,16), (12,15), (11,14)].into_iter().collect();
        let res3: HashSet<(usize, usize)> = vec![(13usize, 15usize), (13,16), (12,15), (14,16), (14,17), (14,18)].into_iter().collect();
        assert_eq!(cuts, vec![res3,res2.clone(),res.clone()]);
        let mut r_max = g.find_r_set_max_with_k(&x, &y, 5).unwrap();
        let cuts = g.enumerate_important_cuts(&mut r_max, &y, &mut HashSet::new(), 5);
        assert_eq!(cuts, vec![res2.clone(),res.clone()]);
    }

    #[test]
    fn skew_edge_multicut_test() {
        let gr = Cursor::new("18 24 0\n2 4 5\n1 10 11\n12 13 14 15\n6\n6\n\
                            7 8 9\n16\n17\n17\n18\n18\n1\n1\n1\n1\n\n\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let ordered_node_pairs: Vec<(usize, usize)> = vec![(0usize, 15usize), (1,16), (2,17)];
        let res: HashSet<(usize, usize)> = vec![(0usize, 1usize), (3,5), (4,5)].into_iter().collect();
        let mut skew_should = SkewEdgeMulticutInstance::new(g.clone(), ordered_node_pairs.clone(), 3);
        let mut skew_not = SkewEdgeMulticutInstance::new(g, ordered_node_pairs, 2);
        let cut = skew_should.solve();
        let no_cut = skew_not.solve();
        assert!(no_cut.is_none());
        assert_eq!(cut, Some(res));
    }

    #[test]
    fn adv_skew_edge_multicut_test() {
        let sk = Cursor::new("21 27 0\n10\n11\n12\n13\n14\n15\n16\n\
                                17\n18\n3\n19\n2 5\n3\n20\n5 7 9\n9\n6\n21\n\
                                1 4\n4 6\n8 6\n");
        let sk_c = Digraph::read_graph(sk);
        assert!(sk_c.is_ok());
        let sk_c = sk_c.unwrap();
        let ordered_node_pairs: Vec<(usize, usize)> = vec![(18usize, 10usize), (19,13), (20,17)];
        let res: HashSet<(usize, usize)> = vec![(2usize, 11usize), (5,14)].into_iter().collect();
        let mut skew = SkewEdgeMulticutInstance::new(sk_c.clone(), ordered_node_pairs.clone(), 2);
        let cut = skew.solve();
        assert_eq!(cut, Some(res));
    }
}
