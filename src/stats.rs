use crate::dfvs_instance::DFVSInstance;
use std::sync::mpsc::Receiver;
use std::time::Instant;
use crate::cust_errors::ProcessingError;
use crate::reduction_rules::Rule;
use std::sync::mpsc::TryRecvError::Disconnected;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RuleStats {
    pub rule: Rule,
    pub reduced_nodes: u64,
    pub reduced_edges: i64,
    pub time_took: u128,
    pub suc_apps: usize,
}

impl RuleStats {

    pub fn new(rule: Rule) -> Self {
        RuleStats {
            rule,
            reduced_nodes: 0,
            reduced_edges: 0,
            time_took: 0,
            suc_apps: 0,
        }
    }

    pub fn add(&mut self, n: u64, m: i64, time: u128) {
        self.reduced_nodes += n;
        self.reduced_edges += m;
        self.time_took += time;
        if n > 0 || m > 0 {
            self.suc_apps += 1;
        }
    }

    pub fn add_mult(&mut self, n: u64, m: i64, time: u128, amount: usize) {
        self.reduced_nodes += n;
        self.reduced_edges += m;
        self.time_took += time;
        self.suc_apps += amount;
    }

    /// Adds the values of `other` to `self`.
    pub fn add_2(&mut self, other: &Self) {
        assert_eq!(self.rule, other.rule);
        self.reduced_edges += other.reduced_edges;
        self.reduced_nodes += other.reduced_nodes;
        self.time_took += other.time_took;
        self.suc_apps += other.suc_apps;
    }

    /// Adds the values of `other` to `self`.
    pub fn add_2_new(&self, other: &Self) -> Self {
        assert_eq!(self.rule, other.rule);
        RuleStats {
            rule: self.rule,
            reduced_edges: self.reduced_edges + other.reduced_edges,
            reduced_nodes: self.reduced_nodes + other.reduced_nodes,
            time_took: self.time_took + other.time_took,
            suc_apps: self.suc_apps + other.suc_apps,
        }
    }

    /// TODO first find out what we load into this.
    /// Adds the values of each element of `other` to the matching elements of `this`. 
    /// This function won't work properly if `this` and `other` dont hold matching elements in the
    /// same order.
    /// 
    /// # Arguments
    ///
    /// * `order` - If set to `true` elements in `one` will be merged first, otherwise `other` will
    /// be preferenced.
    ///
    /// # Panics
    /// Panics if `other` is shorter than `this`.
    pub fn merge_vecs(one: &Vec<Self>, other: &Vec<Self>, order: bool) -> Vec<Self> {
        let mut out = Vec::new();
        let mut e = 0;
        let mut r = 0;
        while e < one.len() && r < other.len() {
            if one[e].rule == other[r].rule {
                out.push(one[e].add_2_new(&other[r]));
                e+=1;
                r+=1;
            } else {
                if order {
                    out.push(one[e]);
                    e+=1;
                } else {
                    out.push(other[r]);
                    r+=1;
                }
            }
        }
        while e < one.len() {
            out.push(one[e]);
            e+=1;
        }
        while r < other.len() {
            out.push(other[r]);
            r+=1;
        }
        return out
    }

}

impl DFVSInstance {

    /// Applies the global lossy2 rule once on the instance and records the running time and the
    /// kill count.
    pub fn apply_global_lossy2_once(&mut self, param: usize) -> RuleStats {
        let mut rs = RuleStats::new(Rule::GlobalLossy2(param));
        let nodes_before: u64 = self.graph.num_nodes() as u64;
        let edges_before: i64 = self.graph.num_edges() as i64;
        let start_time = Instant::now();
        let amount = self.apply_lossy2_global_rule(param);
        rs.add_mult(
            nodes_before - self.graph.num_nodes() as u64,
            edges_before - self.graph.num_edges() as i64,
            start_time.elapsed().as_millis(),
            amount
        );
        rs
    }


    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// Records the running time and kill count, as well as the number of successfull applications
    /// for each rule.
    /// The rules are applied more fine grained to make it easier to interrupt.
    /// The priority order should roughly be chosen by the time consumption of the respective rules. 
    ///
    /// Simple rules have to be the first rules applied.
    /// # Panics 
    /// Panics if `priority_list` contains a forbidden rule
    pub fn exhaustive_fine_rules_stats(&mut self, priority_list: &Vec<Rule>, rec: &Receiver<u8>) -> Result<Vec<RuleStats>, ProcessingError> {
        let mut rule_stats: Vec<RuleStats> = priority_list.iter().map(|rule| RuleStats::new(*rule)).collect();
        'outer: loop {
            for rule_stat in &mut rule_stats {
                let tr = rec.try_recv();
                match tr {
                    Err(Disconnected) => {
                        eprintln!("interrupted since disco");
                        return Err(ProcessingError::OutOfTime);
                    },
                    Ok(_) => {
                        eprintln!("interrupted since interrupt send");
                        return Err(ProcessingError::OutOfTime);
                    },
                    _ => (),
                }
                match rule_stat.rule {
                    Rule::SimpleRules => {
                        let nodes_before: u64 = self.graph.num_nodes() as u64;
                        let edges_before: i64 = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        self.apply_simple_rules();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        let tr = rec.try_recv();
                        match tr {
                            Err(Disconnected) => {
                                eprintln!("interrupted since disco");
                                return Err(ProcessingError::OutOfTime);
                            },
                            Ok(_) => {
                                eprintln!("interrupted since interrupt send");
                                return Err(ProcessingError::OutOfTime);
                            },
                            _ => (),
                        }
                    },
                    Rule::SCC => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_advanced_scc_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Dome => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let start_time = Instant::now();
                        let mut app = false;
                        let edges: Vec<_> = self.graph.weak_edges().collect();
                        for edge in edges {
                            app = self.single_dome_rule(edge);
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            0,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Clique => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let mut app = false;
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_clique_rule(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Core => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let mut app = false;
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_core_rule(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::LinkNode => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let mut app = false;
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_link_node_rules(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Crown => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_crown_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::TwinNodes => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_twin_nodes_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Dominion => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let mut app = false;
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_dominion_rule(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Petal => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let mut app = false;
                        self.compute_and_set_fast_upper(true);
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_petal_rules(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedPetal => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        self.compute_and_set_fast_upper(true);
                        let mut app = false;
                        let nodes: Vec<_> = self.graph.nodes().collect();
                        for node in nodes {
                            app = self.single_advanced_petal_rules(node);
                            if app {break};
                            let tr = rec.try_recv();
                            match tr {
                                Err(Disconnected) => {
                                    eprintln!("interrupted since disco");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                Ok(_) => {
                                    eprintln!("interrupted since interrupt send");
                                    return Err(ProcessingError::OutOfTime);
                                },
                                _ => (),
                            }
                        }
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Lossy(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_lossy_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::SimpleLossy2(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_simple_lossy2_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedLossy2(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_advanced_lossy2_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    _ => panic!(),
                }
            }
            break
        }
        return Ok(rule_stats)
    }

    /// Applies the different rules in the order of `priority_list` each time a rule reduced the instance the function starts from the top.
    /// Records the running time and kill count, as well as the number of successfull applications
    /// for each rule.
    /// The priority order should roughly be chosen by the time consumption of the respective rules. 
    ///
    /// Simple rules have to be the first rules applied.
    ///
    /// # Panics 
    /// Panics if `priority_list` contains a forbidden rule
    pub fn exhaustive_rules_stats(&mut self, priority_list: &Vec<Rule>, rec: Receiver<u8>) -> Result<Vec<RuleStats>, ProcessingError> {
        let mut rule_stats: Vec<RuleStats> = priority_list.iter().map(|rule| RuleStats::new(*rule)).collect();
        'outer: loop {
            for rule_stat in &mut rule_stats {
                let tr = rec.try_recv();
                match tr {
                    Err(Disconnected) => {
                        eprintln!("interrupted since disco");
                        return Err(ProcessingError::OutOfTime);
                    },
                    Ok(_) => {
                        eprintln!("interrupted since interrupt send");
                        return Err(ProcessingError::OutOfTime);
                    },
                    _ => (),
                }
                match rule_stat.rule {
                    Rule::SimpleRules => {
                        let nodes_before: u64 = self.graph.num_nodes() as u64;
                        let edges_before: i64 = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        self.apply_simple_rules();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                    },
                    Rule::SCC => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_advanced_scc_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Dome => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let start_time = Instant::now();
                        let app = self.apply_dome_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            0,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Clique => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_exhaustive_clique_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Core => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_exhaustive_core_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::LinkNode => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_link_node_rules();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Crown => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_crown_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::TwinNodes => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_twin_nodes_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Dominion => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_dominion_rule();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Petal => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        self.compute_and_set_fast_upper(true);
                        let app = self.apply_petal_rules();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedPetal => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        self.compute_and_set_fast_upper(true);
                        let app = self.apply_advanced_petal_rules();
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::Lossy(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_lossy_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::SimpleLossy2(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_simple_lossy2_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    Rule::AdvancedLossy2(q) => {
                        let nodes_before = self.graph.num_nodes() as u64;
                        let edges_before = self.graph.num_edges() as i64;
                        let start_time = Instant::now();
                        let app = self.apply_advanced_lossy2_rules(q);
                        rule_stat.add(
                            nodes_before - self.graph.num_nodes() as u64,
                            edges_before - self.graph.num_edges() as i64,
                            start_time.elapsed().as_millis()
                        );
                        if app {
                            continue 'outer
                        }
                    },
                    _ => panic!(),
                }
            }
            break
        }
        return Ok(rule_stats)
    }

}

