use crate::dfvs_instance::DFVSInstance;
use fxhash::FxHashSet;
use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rule {
    SimpleRules,
    Dome,
    Clique,
    Core,
    SCC,
    LinkNode,
    Crown,
    TwinNodes,
    Dominion,
    Petal,
    AdvancedPetal,
    Lossy(usize),
    SimpleLossy2(usize),
    AdvancedLossy2(usize),
    GlobalLossy2(usize),
}

impl DFVSInstance {

    /// Applies some simple reduction rules. The rules are applied for each node and in no specific
    /// order until no more rules can be applied.
    /// Returns true if at least one reduction has been applied.
    ///
    /// The rules are: 
    /// Rule 0 (implicit): Remove merge multi edges to one. 
    /// Rule 1: Remove loop nodes and add them to the solution. 
    /// Rule 2: Remove sources and sinks.
    /// Rule 3.1: Replace node v with one incoming neighbor u and one outgoing neighbor w by
    /// the edge (u,w).
    /// Rule 3.2: Merge node v with only one incoming neighbor u, back into u.
    /// Rule 3.3: Merge node v with only one outgoing neighbor w, into w.
    pub fn apply_simple_rules(&mut self) -> bool {
        let mut changed = true;
        let mut rounds = 0;
        while changed {
            rounds +=1;
            changed = false;
            let nodes = self.graph.nodes().collect::<Vec<usize>>();
            for node in nodes {
                // Rule 1: Add nodes with loops to the solution, remove it:
                if self.graph.in_neighbors(node).as_ref().expect("`node` is in graph.nodes()").contains(&node) {
                    self.add_to_solution(node).expect("`node` exists");
                    changed = true;
                // Rule 2: Remove sources and sinks:
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 0 || self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 0 {
                    self.remove_node(node).expect("`node` exists");
                    changed = true;
                // Rule 3.1: Replace node v with one incoming neighbor u and one outgoing neighbor w by
                // the edge (u,w):
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 1 && self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let src = self.graph.in_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("Indegree of `node` > 0");
                    let trg = self.graph.out_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("Outdegree of `node` > 0");
                    self.replace(node, (src, trg)).expect("`node` exists");
                    changed = true;
                // Rule 3.2: Merge node v with only one incoming neighbor u, back into u.
                } else if self.graph.in_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let into = self.graph.in_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("`node` has indegree > 0");
                    self.merge_into_back(node, into).expect("Nodes are suited for this operation");
                    changed = true;
                // Rule 3.3: Merge node v with only one outgoing neighbor w, into w.
                } else if self.graph.out_degree(node).expect("`node` is in graph.nodes()") == 1 {
                    let into = self.graph.out_neighbors(node).clone().expect("`node` is in graph.nodes()").into_iter().next().expect("`node` has outdegree > 0");
                    self.merge_into_front(node, into).expect("Nodes are suited for this operation");
                    changed = true;
                }
            }
            // Rule 0: Remove duplicate edges is implicit. 
        }
        rounds>1
    }

    /// Applies the `LinkNode`-rule that contracts nodes of a strict strong degree of 2 (and no
    /// other incident edges) and merges its neighbors `neighs[0]` and `neighs[1]`
    /// (or adds both of them to the solution if they are strongly connected. 
    /// If there exists a strictly weak path between `neighs[0]` and `neighs[1]`, without there
    /// being an edge connecting those two, this rule is not applicable.
    /// "Small clique rule" only needs to be applied if clique rule wasn't applied before this
    /// rules.
    pub fn apply_link_node_rules(&mut self) -> bool {
        let mut changed = false;
        'outer: loop {
            let nodes = self.graph.nodes().collect::<Vec<usize>>();
            for node in nodes {
                if let Some(neighs) = self.graph.is_link_node(node) {
                    if self.graph.strongly_connected(neighs[0], neighs[1]) {
                        // Clique rule has not been applied yet. 
                        self.add_to_solution(neighs[0]).expect("Given node exists");
                        self.add_to_solution(neighs[1]).expect("Given node exists");
                        changed = true;
                        continue 'outer
                    } else {
                        // if there does not exist a edge from neighs[0] to neighs[1]:
                        if !self.graph.in_neighbors(neighs[1]).as_ref().expect("`neighs[1]` exists")
                            .contains(&neighs[0]) {
                            //  check if a weak path from neighs[1] to neighs[0] exist. 
                            let mut weak_in_0 = self.graph.weak_in_neighbors(neighs[0]).expect("`neighs[0]` exists");
                            weak_in_0.remove(&neighs[1]);
                            let mut weak_out_1 = self.graph.weak_out_neighbors(neighs[1]).expect("`neighs[1]` exists");
                            weak_out_1.remove(&neighs[0]);
                            if !(weak_in_0.is_empty() || weak_out_1.is_empty()) {
                                if self.graph.weak_path_exists_between(&weak_out_1, &weak_in_0) {
                                    continue
                                }
                            }
                        }
                        // if there does not exist a edge from neighs[1] to neighs[0]:
                        if !self.graph.in_neighbors(neighs[0]).as_ref().expect("`neighs[0]` exists")
                            .contains(&neighs[1]) {
                            //  check if a weak path from neighs[1] to neighs[0] exist. 
                            let mut weak_in_1 = self.graph.weak_in_neighbors(neighs[1]).expect("`neighs[1]` exists");
                            weak_in_1.remove(&neighs[0]);
                            let mut weak_out_0 = self.graph.weak_out_neighbors(neighs[0]).expect("`neighs[0]` exists");
                            weak_out_0.remove(&neighs[1]);
                            if !(weak_in_1.is_empty() || weak_out_0.is_empty()) {
                                if self.graph.weak_path_exists_between(&weak_out_0, &weak_in_1) {
                                    continue
                                }
                            }
                        }
                        self.contract_link_node(node, &neighs).expect("Nodes are suited for this operation");
                        changed = true;
                        continue 'outer
                    }
                }
            }
            break 'outer
        }
        changed
    }

    /// Applies the `LinkNode`-rule on `node`. If `node` has a strict strong degree of 2 (and no
    /// other incident edges) it is contracted and its neighbors `neighs[0]` and `neighs[1]` are
    /// merged or added to the solution if they are strongly connected. 
    /// If there exists a strictly weak path between `neighs[0]` and `neighs[1]`, without there
    /// being an edge connecting those two, this rule is not applicable.
    /// "Small clique rule" only needs to be applied if clique rule wasn't applied before this
    /// rules.
    ///
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_link_node_rules(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        if let Some(neighs) = self.graph.is_link_node(node) {
            if self.graph.strongly_connected(neighs[0], neighs[1]) {
                // Clique rule has not been applied yet. 
                self.add_to_solution(neighs[0]).expect("Given node exists");
                self.add_to_solution(neighs[1]).expect("Given node exists");
                return true;
            } else {
                // if there does not exist a edge from neighs[0] to neighs[1]:
                if !self.graph.in_neighbors(neighs[1]).as_ref().expect("`neighs[1]` exists")
                    .contains(&neighs[0]) {
                    //  check if a weak path from neighs[1] to neighs[0] exist. 
                    let mut weak_in_0 = self.graph.weak_in_neighbors(neighs[0]).expect("`neighs[0]` exists");
                    weak_in_0.remove(&neighs[1]);
                    let mut weak_out_1 = self.graph.weak_out_neighbors(neighs[1]).expect("`neighs[1]` exists");
                    weak_out_1.remove(&neighs[0]);
                    if !(weak_in_0.is_empty() || weak_out_1.is_empty()) {
                        if self.graph.weak_path_exists_between(&weak_out_1, &weak_in_0) {
                            return false
                        }
                    }
                }
                // if there does not exist a edge from neighs[1] to neighs[0]:
                if !self.graph.in_neighbors(neighs[0]).as_ref().expect("`neighs[0]` exists")
                    .contains(&neighs[1]) {
                    //  check if a weak path from neighs[1] to neighs[0] exist. 
                    let mut weak_in_1 = self.graph.weak_in_neighbors(neighs[1]).expect("`neighs[1]` exists");
                    weak_in_1.remove(&neighs[0]);
                    let mut weak_out_0 = self.graph.weak_out_neighbors(neighs[0]).expect("`neighs[0]` exists");
                    weak_out_0.remove(&neighs[1]);
                    if !(weak_in_1.is_empty() || weak_out_0.is_empty()) {
                        if self.graph.weak_path_exists_between(&weak_out_0, &weak_in_1) {
                            return false
                        }
                    }
                }
                self.contract_link_node(node, &neighs).expect("Nodes are suited for this operation");
                return true
            }
        }
        return false
    }

    /// Traverses over all nodes with a strict strong degree of 3 (and no weak degree) and stores
    /// the neighbors `neighs` of each such nodes either in `connects` if any pair in `neighs` is
    /// connected by a strong edge, or in `un_connects` otherwise. 
    /// If for any set of `neighs` there already exists an entry in `connects`, adds all nodes in
    /// `neighs` to the solution and returns. If `neighs` are already in `un_connects`, stores
    /// `neighs` also in `connects`.
    ///
    /// TODO merge if no edge exists instead of using `un_connects`.
    pub fn apply_twin_nodes_rule(&mut self) -> bool {
        let mut connects = HashSet::new();
        let mut un_connects = HashSet::new();
        let nodes = self.graph.nodes().collect::<Vec<usize>>();
        for node in nodes {
            if let Some(neighs) = self.graph.is_trip_node(node) {
                if connects.contains(&neighs) {
                    self.add_to_solution(neighs[0]).expect("Node exists");
                    self.add_to_solution(neighs[1]).expect("Node exists");
                    self.add_to_solution(neighs[2]).expect("Node exists");
                    return true
                }
                else if un_connects.contains(&neighs) {
                    connects.insert(neighs);
                } else {
                    if self.graph.has_strong_edge(&neighs.iter().copied().collect()) {
                        connects.insert(neighs);
                    } else {
                        un_connects.insert(neighs);
                    }
                }
            }
        }
        false
    }

    /// Finds all strongly connected components and removes the edges between them and single node
    /// components.
    ///
    /// Returns true if at least one edge or node has been removed.
    pub fn apply_scc_rule(&mut self) -> bool {
        let mut change = false;
        let sccs = self.graph.find_strongly_connected_components_iter();
        let matter_sccs: Vec<FxHashSet<usize>> = sccs.into_iter()
            .filter(|scc| {
                if scc.len() == 1 {
                    let elem = scc.iter().next().expect("`scc` holds one element");
                    self.remove_node(*elem).expect("Node exists");
                    change = true;
                    false
                } else {
                    true
                }
            }).collect();
        for i in 0..matter_sccs.len() {
            for j in 0..matter_sccs.len(){
                if i == j {
                    continue
                }
                let edges_from_to: Vec<_> = self.graph.edges_from_to(&matter_sccs[i], &matter_sccs[j]);
                if !edges_from_to.is_empty() {
                    self.remove_edges(edges_from_to).expect("Edges exists");
                    change = true;
                }
            }
        }
        change
    }

    /// Finds all strongly connected components ignoring strong edges and removes the edges
    /// between.
    ///
    /// Returns true if at least one edge has been removed.
    pub fn apply_advanced_scc_rule(&mut self) -> bool {
        let sccs = self.graph.find_weak_strongly_connected_components_iter();
        let mut change = false;
        for i in 0..sccs.len() {
            for j in 0..sccs.len(){
                if i == j {
                    continue
                }
                let edges_from_to: Vec<_> = self.graph.weak_edges_from_to(&sccs[i], &sccs[j]);
                if !edges_from_to.is_empty() {
                    self.remove_edges(edges_from_to).expect("Edges exists");
                    change = true;
                }
            }
        }
        change
    }

    /// Applies either of two rules to the first node one of them appies to.
    /// 1. Nodes which are adjacent to exactly one node disjunct cycle can be contracted.
    /// 2. Nodes which are adjacent to at least `self.upper_bound` + 1 - the lower bound of the
    ///    left over graph many node disjunct cycles
    ///    can be added to the solution.
    pub fn apply_advanced_petal_rules(&mut self) -> bool {
        for node in self.graph.nodes().collect::<Vec<_>>() {
            let (num_petals, left_over) = self.graph.count_petals_left_over(node);
            if num_petals == 1 {
                self.contract_node(node).expect("`node` exists"); 
                return true
            }
            if let Some(upper) = self.effective_upper_bound() {
                let lo_ins = DFVSInstance::new(left_over, None, None);
                let lower = lo_ins.lower_bound_clique_heuristic(&vec![Rule::SimpleRules, Rule::SCC], false);
                // compute lower in left over 
                if upper-lower < num_petals {
                    self.add_to_solution(node).expect("`node` exists");
                    return true
                }
            }
        }
        return false;
    }

    /// Applies either of two rules if any applies to `node`.
    /// 1. Nodes which are adjacent to exactly one node disjunct cycle can be contracted.
    /// 2. Nodes which are adjacent to at least `self.upper_bound` + 1 - the lower bound of the
    ///    left over graph many node disjunct cycles
    ///    can be added to the solution.
    ///
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_advanced_petal_rules(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        let (num_petals, left_over) = self.graph.count_petals_left_over(node);
        if num_petals == 1 {
            self.contract_node(node).expect("`node` exists"); 
            return true
        }
        if let Some(upper) = self.effective_upper_bound() {
            let lo_ins = DFVSInstance::new(left_over, None, None);
            let lower = lo_ins.lower_bound_clique_heuristic(&vec![Rule::SimpleRules, Rule::SCC], false);
            // compute lower in left over 
            if upper-lower < num_petals {
                self.add_to_solution(node).expect("`node` exists");
                return true
            }
        }
        return false;
    }

    /// Applies either of two rules to the first node one of them appies to.
    /// 1. Nodes which are adjacent to exactly one node disjunct cycle can be contracted.
    /// 2. Nodes which are adjacent to at least `self.upper_bound` + 1 many node disjunct cycles
    ///    can be added to the solution.
    pub fn apply_petal_rules(&mut self) -> bool {
        for node in self.graph.nodes().collect::<Vec<_>>() {
            let num_petals = self.graph.count_petals(node);
            if num_petals == 1 {
                self.contract_node(node).expect("`node` exists"); 
                return true
            }
            if let Some(upper) = self.effective_upper_bound() {
                if upper < num_petals {
                    self.add_to_solution(node).expect("`node` exists");
                    return true
                }
            }
        }
        return false;
    }

    /// Applies either of two rules if any applies to `node`.
    /// 1. Nodes which are adjacent to exactly one node disjunct cycle can be contracted.
    /// 2. Nodes which are adjacent to at least `self.upper_bound` + 1 many node disjunct cycles
    ///    can be added to the solution.
    ///
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_petal_rules(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        let num_petals = self.graph.count_petals(node);
        if num_petals == 1 {
            self.contract_node(node).expect("`node` exists"); 
            return true
        }
        if let Some(upper) = self.effective_upper_bound() {
            if upper < num_petals {
                self.add_to_solution(node).expect("`node` exists");
                return true
            }
        }
        return false;
    }

    /// An (hopefully) efficient implementation of the local k-flower rule. 
    /// 
    /// Considers every daisy `D` with at least effective upper bound - effective lower bound + 2 petals
    /// until some node was added to the solution.
    /// Checks if `D` has at least effective upper bound - effective lower bound of G/D + 1 petals,
    /// if so, adds the daisy core to the solution.
    ///
    /// Returns `true` if some node was added to the solution.
    /// For now also returns `false` if execution was interrupted by `self.interrupter`.
    #[deprecated(since = "1.6.0", note = "Better use `apply_advanced_petal_rules()")]
    pub fn apply_local_k_daisy(&mut self) -> bool {
        if self.upper_bound.is_some() && self.lower_bound.is_some() {
            for node in self.graph.nodes().collect::<Vec<_>>() {
                let daisy_petals = self.graph.strong_neighbors(node).expect("`node` exists");
                if daisy_petals.len() < self.effective_upper_bound().expect("`self.upper_bound` is some") - 
                    self.effective_lower_bound().expect("`self.lower_bound` is some") + 2 {
                        continue;
                }
                let mut left_over = self.graph.clone();
                left_over.remove_node(node);
                left_over.remove_nodes(daisy_petals.clone());
                let mut left_over_ins = DFVSInstance::new(left_over, None, None);
                left_over_ins.compute_and_set_lower(false);
                if daisy_petals.len() > self.effective_upper_bound().expect("`self.upper_bound` is some") - 
                    left_over_ins.effective_lower_bound().expect("`self.lower_bound` is some") {
                    self.add_to_solution(node).expect("`node` exists");
                    return true;
                }
            }
        }
        false 
    }

    /// Note: this does not try different matchings. Best to skip this and use the
    /// `via_vertex_cover` heuristic.
    pub fn apply_crown_rule(&mut self) -> bool {
        // find max matching under nodes that are not connected to weak edges 
        let strong_ones: FxHashSet<usize> = self.graph.only_strong_nodes().collect();
        // Get outsiders 
        let (_, outsiders) = self.graph.strong_max_matching_between(&strong_ones, &None);
        // find max matching between outsiders and N(outsiders) M2
        // Get I0 of unmatched nodes in outsiders 
        let n_out = self.graph.open_out_neighbors_of_set(&outsiders);
        let (m2,mut spikes) = self.graph.strong_max_matching_between(&outsiders, &Some(n_out));
        let mut s = 0;
        let mut head = FxHashSet::default();
        while s != spikes.len() {
        // repeat until spikes does not change 
            s = spikes.len();
            //  let H = N(I) 
            head = self.graph.open_out_neighbors_of_set(&spikes);
            spikes.extend(m2.iter().filter_map(|(a, b)| {
                    if head.contains(a) {
                        Some(*b)
                    } else if head.contains(b) {
                        Some(*a)
                    } else {
                        None
                    }
                }));
        }
        if !spikes.is_empty() {
            // add head to the solution.
            self.add_all_to_solution(&head).expect("All nodes in `head` exist");
            // remove spikes 
            self.remove_nodes(&spikes).expect("All nodes in `spikes` exist");
            return true
        }
        false 
    }

    /// Looks for an unconfined vertex and adds it to the solution if one was found.
    /// Returns `true` if a vertex was added to the solution and `false` otherwise.
    pub fn apply_dominion_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<usize>>();
        for node in nodes {
            let mut set: FxHashSet<usize> = vec![node].into_iter().collect();
            // Only the strong neighbors matter here.
            let mut set_closed_n = self.graph.strong_neighbors(node).expect("`node` exists");
            set_closed_n.insert(node);
            loop {
                let set_closed_n_clone = set_closed_n.clone();
                let opt = set_closed_n_clone.iter().filter_map(|neigh| {
                    if set.contains(neigh) {
                        return None
                    }
                    // TODO This can probably be relaxed:
                    if self.graph.min_weak_degree(*neigh).expect("`neigh` exists") != 0 {
                        return None
                    }
                    let nn = self.graph.strong_neighbors(*neigh).expect("`neigh` exists");
                    if nn.intersection(&set).count() != 1 {
                        return None
                    }
                    Some(nn.difference(&set_closed_n).copied().collect::<FxHashSet<usize>>())
                }).min_by_key(|diff| diff.len());
                if let Some(diff) = opt {
                    if diff.is_empty() {
                        self.add_to_solution(node).expect("`node` exists");
                        return true
                    } else if diff.len() == 1 {
                        let s_prime = diff.into_iter().next().expect("`diff.len()` == 1");
                        set.insert(s_prime);
                        set_closed_n.extend(self.graph.strong_neighbors(s_prime).expect("`s_prime` exists"));
                        set_closed_n.insert(s_prime);
                        continue
                    }
                }
                // TODO: do diamond instead of break 
                break
            }
        }
        false
    }

    /// Checks if `node` is an unconfined vertex and adds it to the solution if it is.
    /// Returns `true` if `node` was added to the solution and `false` otherwise.
    ///
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_dominion_rule(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        let mut set: FxHashSet<usize> = vec![node].into_iter().collect();
        // Only the strong neighbors matter here.
        let mut set_closed_n = self.graph.strong_neighbors(node).expect("`node` exists");
        set_closed_n.insert(node);
        loop {
            let set_closed_n_clone = set_closed_n.clone();
            let opt = set_closed_n_clone.iter().filter_map(|neigh| {
                if set.contains(neigh) {
                    return None
                }
                // TODO This can probably be relaxed:
                if self.graph.min_weak_degree(*neigh).expect("`neigh` exists") != 0 {
                    return None
                }
                let nn = self.graph.strong_neighbors(*neigh).expect("`neigh` exists");
                if nn.intersection(&set).count() != 1 {
                    return None
                }
                Some(nn.difference(&set_closed_n).copied().collect::<FxHashSet<usize>>())
            }).min_by_key(|diff| diff.len());
            if let Some(diff) = opt {
                if diff.is_empty() {
                    self.add_to_solution(node).expect("`node` exists");
                    return true
                } else if diff.len() == 1 {
                    let s_prime = diff.into_iter().next().expect("`diff.len()` == 1");
                    set.insert(s_prime);
                    set_closed_n.extend(self.graph.strong_neighbors(s_prime).expect("`s_prime` exists"));
                    set_closed_n.insert(s_prime);
                    continue
                }
            }
            // TODO: do diamond instead of break 
            break
        }
        false
    }

    /// Applies the clique rule: 
    /// Let `C` be an induced strongly connected clique and `n` a node of `C` with n incoming-
    /// or no outgoing neighbors that are not in `C`. Then all nodes of `C`, that are not `n` can be added
    /// to the solution then all nodes in `C` can be removed.
    ///
    /// Better use `apply_exhaustive_clique_rule()`
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_clique_rule()")]
    pub fn apply_clique_rule(&mut self) -> bool {
        let mut greedy_clique = self.graph.greedy_max_clique();
        let mut found_node = None;
        for node in &greedy_clique {
            // no in or out degree?
            if self.graph.min_direct_degree_outside(*node, &greedy_clique).expect("`node` is part of the cliqu") == 0 {
                found_node = Some(*node);
                break;
            }
        }
        if let Some(node) = found_node {
            greedy_clique.remove(&node);
            self.remove_node(node).expect("`node` exists");
            self.add_all_to_solution(&greedy_clique).expect("All nodes in `greedy_clique` exist");
            return true
        }
        false
    }

    /// Applies the clique rule (see documentation) by first finding a greedy 
    /// clique with `self.graph.greedy_max_clique()`, and then going through all the nodes of the
    /// graph, checking whether a node is found which forms a clique with part if the greedy
    /// clique, while having only incoming or outgoing neighbors outisde that clique: 
    ///
    /// Better use `apply_exhaustive_clique_rule()`
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_clique_rule()")]
    pub fn apply_advanced_clique_rule(&mut self) -> bool {
        let greedy_cluster = self.graph.greedy_max_clique();
        let mut addable_cluster = FxHashSet::default();
        let mut removable_node = None;
        for node in self.graph.nodes() {
            if self.graph.min_direct_degree_outside(node, &greedy_cluster).expect("`node` is in `.nodes()`") == 0 {
                let into_sol = self.graph.strong_neighbors(node).expect("`node` is in `.nodes()`").clone();
                if !into_sol.is_empty() && self.graph.min_direct_degree_outside(node, &into_sol).expect("`node` exists") == 0 {
                    addable_cluster = into_sol;
                    removable_node = Some(node);
                    break;
                }
            }
        }
        if let Some(node) = removable_node {
            self.remove_node(node).expect("`node` exists");
            self.add_all_to_solution(&addable_cluster).expect("All nodes in the set exist");
            return true
        }
        false
    }

    /// Applies the clique rule (see documentation) exhaustively by checking for each node
    /// `v` in `self.graph` if `v` has a minimum weak degree of 0. If so, checks if the
    /// strong neighborhood of `v` is a clique. If this is given as well, adds all strong neighbors
    /// of `v` to the solution and removes the neighborhood as well as `v`.
    /// Returns `true` if anything happened.
    pub fn apply_exhaustive_clique_rule(&mut self) -> bool {
        // Go through all the nodes 
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for node in nodes {
            if self.graph.min_weak_degree(node) != Some(0) {
                continue
            }
            // Get strong neighborhood 
            let strong_neighborhood = self.graph.strong_neighbors(node).expect("`node` exists");
            // Check if strong neighborhood is a cluster 
            if self.graph.is_cluster(&strong_neighborhood) {
                self.remove_node(node).expect("`node` exists");
                self.add_all_to_solution(&strong_neighborhood).expect("All nodes in the set exist");
                return true
            }
        }
        false 
    }

    /// Applies the clique rule (see documentation) by checking for `node` in`self.graph` if `node` has a
    /// minimum weak degree of 0. If so, checks if the strong neighborhood of `node` is a clique. 
    /// If this is given as well, adds all strong neighbors of `node` to the solution and removes the 
    /// neighborhood as well as `node`.
    /// Returns `true` if anything happened.
    ///
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_clique_rule(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        if self.graph.min_weak_degree(node) != Some(0) {
            return false
        }
        // Get strong neighborhood 
        let strong_neighborhood = self.graph.strong_neighbors(node).expect("`node` exists");
        // Check if strong neighborhood is a cluster 
        if self.graph.is_cluster(&strong_neighborhood) {
            self.remove_node(node).expect("`node` exists");
            self.add_all_to_solution(&strong_neighborhood).expect("All nodes in the set exist");
            return true
        }
        false 
    }

    /// Applies the core rule (see documentation) by greedily finding a near maximal clique
    /// `clique`, then looking for the node `v` with the least minimum amount of direct neighbors `
    /// neighs` outside of the clique and then checking, if the `neighs` have a common clique
    /// core in `clique`. This core is then added to the solution and removed from the graph.
    /// Returns `true` if a core was added to the solution.
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_core_rule()")]
    pub fn apply_core_clique_rule(&mut self) -> bool {
        // Use heuristic to find max cluster 
        let mut greedy_cluster = self.graph.greedy_max_clique();
        // Pick `node` with lowest min directed degree (outside)
        if !greedy_cluster.is_empty(){
            // TODO: why do we need the clone here?
            let greedy_clone = greedy_cluster.clone();
            let (min_node, its_min_neighbors) = greedy_clone.iter()
                .map(|node| (node, self.graph.min_direct_neighbors_outside(*node, &greedy_cluster).expect("`node` exists")))
                .min_by_key(|(_, neighs)| neighs.len()).expect("`greedy_cluster` is not empty");
            // Remove `node` from cluster
            greedy_cluster.remove(min_node);
            // Intesect the strong neighbors of all min directed neighbors of `node` with max cluster. 
            for node in its_min_neighbors {
                greedy_cluster = self.graph.strong_neighbors_in(node, &greedy_cluster).expect("`node` exists");
            }
        }
        // Add intersection to solution.
        self.add_all_to_solution(&greedy_cluster).expect("All nodes in the set exist");
        !greedy_cluster.is_empty() 
    }

    /// Applies the core rule (see documentation) on the maximal daisy.
    /// Returns `true` if daisy core was added to the solution.
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_core_rule()")]
    pub fn apply_daisy_core_rule(&mut self) -> bool {
        // Get max daisy
        if let Some((daisy_core, daisy_leaves)) = self.graph.get_max_strong_degree_node_and_neighbors() {
            // Create union
            let mut daisy = daisy_leaves.clone();
            daisy.insert(daisy_core);
            // find node in daisy that has a min_direct_degree_outside union of 0. 
            for node in daisy_leaves {
                if self.graph.min_direct_degree_outside(node, &daisy) == Some(0) {
                    self.add_to_solution(daisy_core).expect("This node exists");
                    return true
                }
            }
        }
        false
    }

    /// Applies the core rule on all possible daisies.
    /// Returns `true` if a core was added to the solution.
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_core_rule()")]
    pub fn apply_exhaustive_daisy_core_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for daisy_core in nodes {
            if let Some(daisy_leaves) = self.graph.strong_neighbors(daisy_core) {
                // Create union
                let mut daisy = daisy_leaves.clone();
                daisy.insert(daisy_core);
                // find node in daisy that has a min_direct_degree_outside union of 0. 
                for node in daisy_leaves {
                    if self.graph.min_direct_degree_outside(node, &daisy) == Some(0) {
                        self.add_to_solution(daisy_core).expect("This node exists");
                        return true
                    }
                }
            }
        }
        false
    }

    /// Applies the core rule on the node with the minimal amount of minimal directed neighbors and
    /// looks for a core under those neighbors.
    /// Returns true if something happened.
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_core_rule()")]
    pub fn apply_min_direct_core_rule(&mut self) -> bool {
        // Get min direct degree node
        if let Some((node, min_neighbors)) = self.graph.get_min_min_direct_degree_node_and_neighbors() {
            // Get strong neighbors of `node`:
            if let Some(strong_neighbors) = self.graph.strong_neighbors(node) {
                // find core of `min_neighbors` in `strong_neighbors`: 
                let mut core = strong_neighbors;
                for nigh in min_neighbors {
                    core = self.graph.strong_closed_neighbors_in(nigh, &core).expect("`node` exists");
                    if core.is_empty() {
                        return false
                    }
                }
                for cs in core {
                    self.add_to_solution(cs).expect("This node exists");
                }
                return true
            }
        }
        false
    }

    /// Applies the core rule for all the nodes in the graph, where we look for a core under the
    /// minimal directed neighbors of the node.
    /// Returns true if something happened.
    #[deprecated(since = "1.6.0", note = "Better use `apply_exhaustive_core_rule()")]
    pub fn apply_exhaustive_min_direct_core_rule(&mut self) -> bool {
        // Get min direct degree node
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        'outer: for node in nodes {
            if let Some(min_neighbors) = self.graph.min_direct_neighbors(node) {
                // Get strong neighbors of `node`:
                if let Some(strong_neighbors) = self.graph.strong_neighbors(node) {
                    // find core of `min_neighbors` in `strong_neighbors`: 
                    let mut core = strong_neighbors.clone();
                    for nigh in min_neighbors {
                        core = self.graph.strong_closed_neighbors_in(nigh, &core).expect("`node` exists");
                        if core.is_empty() {
                            continue 'outer
                        }
                    }
                    for cs in core {
                        self.add_to_solution(cs).expect("This node exists");
                    }
                    return true
                }
            }
        }
        false
    }

    /// Applies the core rule for all the nodes in the graph, where those nodes are the fix point
    /// of the core.
    ///
    /// Returns true if something happened.
    pub fn apply_exhaustive_core_rule(&mut self) -> bool {
        let nodes = self.graph.nodes().collect::<Vec<_>>();
        for node in nodes {
            let strong_neighbors = self.graph.strong_neighbors(node).expect("`node` exists");
            if strong_neighbors.is_empty() {
                continue;
            }
            let in_neighbors = self.graph.in_neighbors(node).as_ref().expect("`node` exists");
            let out_neighbors = self.graph.out_neighbors(node).as_ref().expect("`node` exists");
            let mut core = strong_neighbors.clone();
            for nigh in in_neighbors {
                core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
                if core.is_empty() {
                    break;
                }
            }
            if !core.is_empty() {
                for cn in core {
                    self.add_to_solution(cn).expect("This node exists");
                }
                return true
            }
            let mut core = strong_neighbors.clone();
            for nigh in out_neighbors {
                core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
                if core.is_empty() {
                    break;
                }
            }
            if !core.is_empty() {
                for cn in core {
                    self.add_to_solution(cn).expect("This node exists");
                }
                return true
            }
        }
        false
    }

    /// Applies the core rule for a given node `node` in the graph, where `node` is the fix point
    /// of the core.
    /// Returns true if something happened.
    /// 
    /// # Panics 
    /// Panics if `node` is not in `self.graph`.
    pub fn single_core_rule(&mut self, node: usize) -> bool {
        assert!(self.graph.has_node(node));
        let strong_neighbors = self.graph.strong_neighbors(node).expect("`node` exists");
        if strong_neighbors.is_empty() {
            return false;
        }
        let in_neighbors = self.graph.in_neighbors(node).as_ref().expect("`node` exists");
        let out_neighbors = self.graph.out_neighbors(node).as_ref().expect("`node` exists");
        let mut core = strong_neighbors.clone();
        for nigh in in_neighbors {
            core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
            if core.is_empty() {
                break;
            }
        }
        if !core.is_empty() {
            for cn in core {
                self.add_to_solution(cn).expect("This node exists");
            }
            return true
        }
        let mut core = strong_neighbors.clone();
        for nigh in out_neighbors {
            core = self.graph.strong_closed_neighbors_in(*nigh, &core).expect("`node` exists");
            if core.is_empty() {
                break;
            }
        }
        if !core.is_empty() {
            for cn in core {
                self.add_to_solution(cn).expect("This node exists");
            }
            return true
        }
        false
    }

    /// Applies the dome rule that removes edges that are not part of a minimal cycle.
    /// Returns true if something changed.
    pub fn apply_dome_rule(&mut self) -> bool {
        let weak_edges = self.graph.weak_edges().collect::<Vec<_>>();
        let mut change = false;
        for (src, trg) in weak_edges {
            let weak_ins = self.graph.weak_in_neighbors(src).expect("`src` exists");
            let ins = self.graph.in_neighbors(trg).as_ref().expect("`trg` exists");
            if weak_ins.is_subset(ins) {
                self.remove_edge(&(src, trg)).expect("This edge exists");
                change = true;
                continue
            }
            let weak_outs = self.graph.weak_out_neighbors(trg).expect("`src` exists");
            let outs = self.graph.out_neighbors(src).as_ref().expect("`trg` exists");
            if weak_outs.is_subset(outs) {
                self.remove_edge(&(src, trg)).expect("This edge exists");
                change = true;
            }
        }
        change
    }

    /// Applies the dome rule. Checks if `edge` is not part of a minimal cycle and in that case
    /// removes it.
    /// Returns true if something changed.
    ///
    /// # Panics
    /// Panics if `edge` does not exist in `self.graph`.
    ///
    /// This rule can be applied for all weak_edges in `self.graph` before returning to the simple
    /// rules. TODO: double check
    pub fn single_dome_rule(&mut self, edge: (usize, usize)) -> bool {
        assert!(self.graph.is_weak_edge(edge));
        let (src, trg) = edge;
        let weak_ins = self.graph.weak_in_neighbors(src).expect("`src` exists");
        let ins = self.graph.in_neighbors(trg).as_ref().expect("`trg` exists");
        if weak_ins.is_subset(ins) {
            self.remove_edge(&(src, trg)).expect("This edge exists");
            return true;
        }
        let weak_outs = self.graph.weak_out_neighbors(trg).expect("`src` exists");
        let outs = self.graph.out_neighbors(src).as_ref().expect("`trg` exists");
        if weak_outs.is_subset(outs) {
            self.remove_edge(&(src, trg)).expect("This edge exists");
            return true
        }
        false 
    }

    /// Applies a lossy kernelization rule that adds strong cliques of size > `quality` to the
    /// solution. 
    /// By the repeated application of this rule with parameter `quality`, the size of an optimal solution for the resulting kernel can become at worst 1 + 1/`quality` times as large.
    pub fn apply_lossy_rules(&mut self, quality: usize) -> bool {
        let greed_clique = self.graph.greedy_max_clique();
        let len = greed_clique.len();
        if len > quality {
            self.add_all_to_solution(&greed_clique).expect("all exist");
            return true
        }
        false
    }

    /// Applies a lossy kernelization rule that contracts nodes with at most `max(out_degree, in_degree)` <= `quality`.
    /// By the repeated application of this rule with parameter `quality`, the size of an optimal solution for the resulting kernel can become at worst `quality` times as large.
    ///
    /// Attention: This one does not work!
    pub fn apply_simple_lossy2_rules(&mut self, quality: usize) -> bool {
        if let Some((node,neighs)) = self.graph.get_min_min_direct_degree_node_and_neighbors() {
            if neighs.len() <= quality {
                self.contract_node(node).expect("`node` exists");
                return true
            }
            return false 
        } 
        false
    }

    /// Applies a lossy kernelization rule that contracts nodes incident to at most `quality` many
    /// petals.
    /// By the repeated application of this rule with parameter `quality`, the size of an optimal solution for the resulting kernel can become at worst `quality` times as large.
    ///
    /// Attention: This one does not work!
    pub fn apply_advanced_lossy2_rules(&mut self, quality: usize) -> bool {
        for node in self.graph.nodes().collect::<Vec<_>>() {
            let num_petals = self.graph.count_petals(node);
            if num_petals <= quality {
                self.contract_node(node).expect("`node` exists"); 
                return true
            }
        }
        return false;
    }

    /// Applies a lossy kernelization rule that contracts nodes incident to at most `quality` many
    /// petals.
    /// This rule is applied on multiple different cycle disjunct nodes and with this maintaining a
    /// approximation factor of `quality`. 
    ///
    /// If after the application of this rule another lossy rule is applied, the approximation
    /// factor will not be stable. In the worst case the factor can become `quality` times the
    /// quality of the applied lossy rules. E.g. if after this rule the lossy1 rules are applied
    /// with a parameter of 1 (which normally leads to an approximation factor of 2) the new
    /// approximation factor will be 4 in the worst case.
    pub fn apply_lossy2_global_rule(&mut self, quality: usize) -> bool {
        let mut nodes: FxHashSet<_> = self.graph.nodes().collect();
        let mut something = false;
        while !nodes.is_empty() {
            let node = *nodes.iter().next().expect("not empty");
            nodes.remove(&node);
            let (num_petals, to_remove) = self.graph.count_petals_give_petals(node);
            if num_petals <= quality {
                self.contract_node(node).expect("`node` exists"); 
                nodes = nodes.difference(&to_remove).copied().collect();
                something = true;
            }
            
        }
        return something;
    }

    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// The priority order should roughly be chosen by the time consumption of the respective rules. 
    ///
    /// Simple rules have to be the first rules applied.
    ///
    /// # Panics
    /// Panics if rules are used that are not supposed to be used here. For example the
    /// `GlobalLossy2` rule.
    ///
    /// TODO: short k flower (only on a few nodes with high strong degree).
    pub fn exhaustive_reductions(&mut self, priority_list: &Vec<Rule>) {
        assert_eq!(priority_list[0], Rule::SimpleRules);
        'outer: loop {
            for rule in priority_list {
                match rule {
                    Rule::SimpleRules => {
                        self.apply_simple_rules();
                    },
                    Rule::SCC => {
                        if self.apply_advanced_scc_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Dome => {
                        if self.apply_dome_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Clique => {
                        if self.apply_exhaustive_clique_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Core => {
                        if self.apply_exhaustive_core_rule() {
                            continue 'outer
                        }
                    },
                    Rule::LinkNode => {
                        if self.apply_link_node_rules() {
                            continue 'outer
                        }
                    },
                    Rule::Crown => {
                        if self.apply_crown_rule() {
                            continue 'outer
                        }
                    },
                    Rule::TwinNodes => {
                        if self.apply_twin_nodes_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Dominion => {
                        if self.apply_dominion_rule() {
                            continue 'outer
                        }
                    },
                    Rule::Petal => {
                        self.compute_and_set_fast_upper(true);
                        if self.apply_petal_rules() {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedPetal => {
                        self.compute_and_set_fast_upper(true);
                        if self.apply_advanced_petal_rules() {
                            continue 'outer
                        }
                    },
                    Rule::Lossy(q) => {
                        if self.apply_lossy_rules(*q) {
                            continue 'outer
                        }
                    },
                    Rule::SimpleLossy2(q) => {
                        if self.apply_simple_lossy2_rules(*q) {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedLossy2(q) => {
                        if self.apply_advanced_lossy2_rules(*q) {
                            continue 'outer
                        }
                    },
                    _ => panic!("Other rules should not be used here!"),
                }
            }
            break
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::digraph::Digraph;
    use crate::dfvs_instance::DFVSInstance;

    #[test]
    fn reduction_test() {
        let gr = Cursor::new("8 12 0\n2 3\n\n2 5\n1\n2 7\n1 4 7\n8 1\n\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let graph = g.unwrap();
        let mut instance = DFVSInstance::new(graph, None, None);
        instance.apply_simple_rules();
        assert_eq!(instance.graph.num_nodes(), 0);
        assert_eq!(instance.solution, vec![0usize].into_iter().collect::<FxHashSet<usize>>());
    }


    #[test]
    fn scc_reduction_test() {
        let gr = Cursor::new("11 19 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7 8 9\n8 9\n9\n11\n8\n10\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut g_check = g.clone();
        let mut instance = DFVSInstance::new(g, None, None);
        instance.apply_scc_rule();
        g_check.remove_node(0);
        g_check.remove_node(5);
        g_check.remove_node(6);
        assert_eq!(instance.graph, g_check);
        let gr = Cursor::new("9 14 0\n2\n3\n4 5\n2 6 7\n4 6 7\n7\n8\n9\n6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut g_check = g.clone();
        let mut instance = DFVSInstance::new(g, None, None);
        instance.apply_scc_rule();
        g_check.remove_node(0);
        g_check.remove_edge(&(3,5));
        g_check.remove_edge(&(3,6));
        g_check.remove_edge(&(4,5));
        g_check.remove_edge(&(4,6));
        let instance_check = DFVSInstance::new(g_check, None, None);
        assert_eq!(instance.graph, instance_check.graph);
    }


    // #[test]
    // fn clique_rule_test() {
    //     let gr = Cursor::new("7 22 0\n2\n1 3 4 5\n2 4 5\n2 3 5 6 7\n2 3 4 7\n4 7\n4 6 5\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(instance.apply_clique_rule());
    //     assert_eq!(instance.solution, vec![1,3,4].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 3);
    // }

    // #[test]
    // fn adv_clique_rule_test() {
    //     let gr = Cursor::new("7 26 0\n2 3 4 6\n1 3 4 7\n1 2 4 5 7\n1 2 3 5 6\n3 4\n1 4 5 7\n2 3 5 6\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(instance.apply_advanced_clique_rule());
    //     assert_eq!(instance.solution, vec![2,3].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 4);
    // }

    #[test]
    fn exh_clique_rule_test() {
        let gr = Cursor::new("7 27 0\n2 3 4 5 6\n1 3 4 7\n1 2 4 7\n1 2 3 6\n6 7\n1 4 5 7\n2 3 5 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_clique_rule());
        assert_eq!(instance.solution, vec![5,6].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);
    }

    #[test]
    fn adv_scc_rule_test() {
        let gr = Cursor::new("6 13 0\n2 3\n3 5\n1 2 6\n3 5\n4 6\n4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_advanced_scc_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 11);
    }

    #[test]
    fn adv_scc_plus_dome_rule_test() {
        let gr = Cursor::new("6 13 0\n3 5\n1 3\n1 4 6\n2 5\n2 6\n4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(!instance.apply_advanced_scc_rule());
        assert!(instance.apply_dome_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 10);
    }

    // #[test]
    // fn core_clique_rule_test() {
    //     let gr = Cursor::new("8 21 0\n2 3 8\n1 3 5 6\n1 2 7 8\n1\n1 2\n2 7\n3 6\n3\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(instance.apply_core_clique_rule());
    //     assert_eq!(instance.solution, vec![2].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 7);
    // }

    // #[test]
    // fn daisy_core_rule_test() {
    //     let gr = Cursor::new("5 13 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(instance.apply_daisy_core_rule());
    //     assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 4);
    // }

    // #[test]
    // fn exh_daisy_core_rule_test() {
    //     let gr = Cursor::new("11 29 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n8 9 10 11\n8 9 10 11\n6 7\n6 7\n6 7\n6 7\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(!instance.apply_daisy_core_rule());
    //     assert!(instance.apply_exhaustive_daisy_core_rule());
    //     assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 10);
    // }

    // #[test]
    // fn min_direct_core_rule_test() {
    //     let gr = Cursor::new("5 13 0\n2 4 5\n1 3 5\n1 2\n1 3 5\n2 3 4\n");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(instance.apply_min_direct_core_rule());
    //     assert_eq!(instance.solution, vec![1].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 4);
    // }

    // #[test]
    // fn exh_min_direct_core_rule_test() {
    //     let gr = Cursor::new("8 25 0\n2 3 4 5\n1 3 4 5\n1 2 4 5\n2 3 6\n2 3 8\n1 7\n1 6 8\n1 7");
    //     let g = Digraph::read_graph(gr);
    //     assert!(g.is_ok());
    //     let g = g.unwrap();
    //     let mut instance = DFVSInstance::new(g, None, None);
    //     assert!(!instance.apply_min_direct_core_rule());
    //     assert!(instance.apply_exhaustive_min_direct_core_rule());
    //     assert_eq!(instance.solution, vec![1, 2].into_iter().collect::<FxHashSet<usize>>());
    //     // remaining graph size
    //     assert_eq!(instance.graph.num_nodes(), 6);
    // }

    #[test]
    fn exhaustive_core_rule_test() {
        let gr = Cursor::new("8 21 0\n2 3 8\n1 3 5 6\n1 2 7 8\n1\n1 2\n2 7\n3 6\n3\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 7);

        let gr = Cursor::new("5 13 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);

        let gr = Cursor::new("11 29 0\n2 3 4\n1 5\n1 2 4\n1 5\n2 3 4\n8 9 10 11\n8 9 10 11\n6 7\n6 7\n6 7\n6 7\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![0].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 10);

        let gr = Cursor::new("5 13 0\n2 4 5\n1 3 5\n1 2\n1 3 5\n2 3 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![1].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 4);

        let gr = Cursor::new("8 25 0\n2 3 4 5\n1 3 4 5\n1 2 4 5\n2 3 6\n2 3 8\n1 7\n1 6 8\n1 7");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_exhaustive_core_rule());
        assert_eq!(instance.solution, vec![1, 2].into_iter().collect::<FxHashSet<usize>>());
        // remaining graph size
        assert_eq!(instance.graph.num_nodes(), 6);
    }

    #[test]
    fn only_adv_scc_rule_test() {
        let gr = Cursor::new("8 17 0\n2 6\n3 7\n1 4 8\n7 5\n8 6\n4 1\n2 4\n3 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(!instance.apply_dome_rule());
        assert!(instance.apply_advanced_scc_rule());
        // remaining graph size
        assert_eq!(instance.graph.num_edges(), 16);
    }

    #[test]
    fn local_k_daisy_rule_test() {
        let gr = Cursor::new("6 8 0\n2 3 4\n1\n1\n1\n6\n5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, Some(3), Some(2));
        assert!(instance.apply_advanced_petal_rules());
        assert_eq!(instance.solution, vec![0].into_iter().collect());
    }

    #[test]
    fn link_node_rule_test() {
        let gr = Cursor::new("7 15 0\n2 7\n1 3\n2 4\n1 3 5\n4 6\n\
                             5 7\n1 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_link_node_rules());
        assert!(instance.rebuild_section().is_ok());
        assert!(instance.apply_link_node_rules());
        instance.finallize_solution(); 
        assert_eq!(instance.solution.len(), 4);
    }

    #[test]
    fn crown_rule_test() {
        //let gr = Cursor::new("12 30 0\n4 5\n5 6\n6 7\n1 5 8 12\n1 2 4 6\n\
        //                     2 3 5 7\n3 6 9 11 12\n4 7 9\n7 8\n8\n10\n4 7\n");
        let gr = Cursor::new("6 9 0\n4\n4\n4\n1 2 3 6\n4\n5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_crown_rule());
        assert_eq!(instance.graph.num_nodes(), 2);
        assert_eq!(instance.solution.len(), 1);
    }

    #[test]
    fn advances_link_node_test() {
        let gr = Cursor::new("13 20 0\n2 3\n1\n1 6 7 8\n2\n2\n7\n\
                              12\n10\n5 10\n9\n4 13\n11\n12\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        let src = vec![5, 6, 7].into_iter().collect();
        let trg = vec![3, 4].into_iter().collect();
        assert!(instance.graph.weak_path_exists_between(&src, &trg));
        assert!(!instance.apply_link_node_rules());
        instance.graph.add_edge((1, 2));
        assert!(instance.apply_link_node_rules());
        instance.apply_simple_rules();
        instance.finallize_solution(); 
        assert_eq!(instance.solution.len(), 3);
        assert!(instance.solution.contains(&0));
    }

    #[test]
    fn twin_node_rule_test() {
        let gr = Cursor::new("5 15 0\n3 4 5\n3 4 5\n1 2 5\n1 2 3\n1 2 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(!instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 0);
        let gr = Cursor::new("6 18 0\n3 4 5\n3 4 5\n1 2 6\n1 2 6\n1 2 6\n\
                              3 4 5\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 3);
        let gr = Cursor::new("5 14 0\n3 4 5\n3 4 5\n1 2\n1 2 5\n1 2 4\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_twin_nodes_rule());
        assert_eq!(instance.solution.len(), 3);
    }

    #[test]
    fn dominion_rule_test() {
        let gr = Cursor::new("14 47 0\n5 6 7\n3 4 7\n2 5 8\n2 5 7\n1 3 4 6 14\n\
                             1 5 7 14\n1 2 4 6 11\n3 9 11 12 13\n10\n8\n7 8 12 13\n\
                             8 11 14\n8 11 14\n5 6 12 13\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
        assert!(instance.apply_dominion_rule());
    }

    #[test]
    fn petal_test() {
        let gr = Cursor::new("10 22 0\n3 5\n4 8\n2 5\n1 8\n7 9\n\
                             1 4 10\n2 3 9\n6 10\n5 7\n6 8\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        instance.compute_and_set_fast_upper(false);
        assert!(instance.apply_petal_rules());
        assert!(!instance.apply_petal_rules());
        let gr = Cursor::new("9 22 0\n3 5\n4 8\n2 5\n1 8\n7 9\n\
                             1 4 9\n2 3 9\n6 9\n5 7 6 8\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        instance.compute_and_set_fast_upper(false);
        assert!(instance.apply_advanced_petal_rules());
        assert!(instance.apply_advanced_petal_rules());
        assert!(instance.apply_advanced_petal_rules());
        assert!(instance.apply_advanced_petal_rules());
        assert_eq!(instance.graph.num_nodes(), 5);
    }

    #[test]
    fn lossy_test() {
        let gr = Cursor::new("10 34 0\n2 3 7 9\n1 3 4 5\n1 2 4 5\n\
                             2 3 5 8 10\n2 3 4 6\n5 7 8\n1 6 8\n\
                             4 6 7\n1 10\n4 9\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let mut instance = DFVSInstance::new(g, None, None);
        let mut clone = instance.clone();
        let mut clone2 = instance.clone();
        let mut moar = 0;
        while instance.apply_lossy_rules(1) {
            moar += 1;
        }
        assert_eq!(moar, 3);
        assert_eq!(instance.graph.num_edges(), 0);
        assert_eq!(instance.solution.len(), 6+3);
        let mut moar = 0;
        while clone.apply_lossy_rules(2) {
            moar += 1;
        }
        assert_eq!(moar, 2);
        clone.apply_simple_rules();
        assert_eq!(clone.graph.num_nodes(), 0);
        assert_eq!(clone.solution.len(), 6+2);
        let mut moar = 0;
        while clone2.apply_lossy_rules(3) {
            moar += 1;
        }
        assert_eq!(moar, 1);
        clone2.apply_simple_rules();
        assert_eq!(clone2.graph.num_nodes(), 0);
        assert_eq!(clone2.solution.len(), 6+1);
    }

}
