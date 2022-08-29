//!
//! This binary is only meant to run experiments and record their results.
//! For each input graph, several output files are created:
//! 1. A csv with general statistics of the graph. Such as number of nodes, number of edges, some
//! upper bound and some lower bound.
//! 2. Two csv's containing statistics of most of the exact rules implemented and the resulting kernels. 
//! One with the `AdvancedPetal`-Rule and one with the `QuickAdvancedPetal`-Rule.
//! 3. Two csv's containing statistics of the essential and the lossy1 rule and the resulting kernel. Again with the two versions of the petal rules.
//! 4. Two csv's containing statistics of the essential, the lossy1 rule and a single application 
//! of the lossy2 rule. After the first (and only) application of the `GlobalLossy2`-Rule 
//! most of the exact rules are exhaustivly applied, and the resulting kernel.
//! 5. Two csv's containing all the rules of 4. with the rules of 3. afterwards.
//!

use std::error;
use std::time::{Duration,Instant};
use std::sync::mpsc::channel;
use std::sync::mpsc::{SendError, Receiver, Sender};
use std::ffi::OsString;
use clap::{Arg, Command};
use std::path::PathBuf;
use std::fs::File;
use std::thread;
use std::thread::JoinHandle;
use std::io::{BufReader, Write};
use std::fmt::Display;

use dfvs_solver::{digraph::Digraph,  dfvs_instance::DFVSInstance, reduction_rules::Rule, stats::RuleStats};

/// Custom Error for the threads.
#[derive(Debug)]
enum ThreadErr {
    SendError(SendError<i8>),
    IoError(std::io::Error),
}

impl Display for ThreadErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThreadErr::SendError(send_error) => 
                write!(f, "{}", send_error),
            ThreadErr::IoError(io_error) => 
                write!(f, "{}", io_error),
        }
    }
}

impl std::error::Error for ThreadErr {}

impl From<SendError<i8>> for ThreadErr {
    fn from(err: SendError<i8>) -> Self {
        ThreadErr::SendError(err)
    }
}

impl From<std::io::Error> for ThreadErr {
    fn from(err: std::io::Error) -> Self {
        ThreadErr::IoError(err)
    }
}

pub fn main() -> Result<(), Box<dyn error::Error>> {

    // CLI stuff
    let m = Command::new("experiments")
        .arg(Arg::new("files")
             .takes_value(true)
             .multiple_values(true)
             .short('f'))
        .arg(Arg::new("dest")
             .required(true)
             .takes_value(true)
             .short('d'))
        .get_matches();
    // Get, as input all public instances. 
    let files: Vec<PathBuf> = m.values_of("files").unwrap().map(|p| PathBuf::from(p)).collect();
    let ultimate = Duration::from_secs(3600);
    let dest: &str = m.value_of("dest").unwrap();

    // Rule priority sets
    let priorities_org = vec![
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::QuickAdvancedPetal],
        vec![Rule::SimpleRules, Rule::LossyClique(1), Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::LossyClique(1), Rule::Dome, Rule::SCC, Rule::QuickAdvancedPetal],
    ];

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/init_stats.csv",dest))?,
        File::create(format!("{}/kern_rules.csv",dest))?,
        File::create(format!("{}/kern_rules_qp.csv",dest))?,
        File::create(format!("{}/sim_rules_l1.csv",dest))?,
        File::create(format!("{}/sim_rules_l1_qp.csv",dest))?,
        File::create(format!("{}/sim_rules_l_all_rules.csv",dest))?,
        File::create(format!("{}/sim_rules_l_all_rules_qp.csv",dest))?,
        File::create(format!("{}/sim_rules_l4_all_rules.csv",dest))?,
        File::create(format!("{}/sim_rules_l4_all_rules_qp.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, n, m, upper_bound, t_upper, lower_bound, t_lower")?;
    writeln!(&mut out_files[1], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[2], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_qap, n_qap, m_qap")?;
    writeln!(&mut out_files[3], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[4], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_qap, n_qap, m_qap")?;
    writeln!(&mut out_files[5], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino")?;
    writeln!(&mut out_files[6], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_qap, n_qap, m_qap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino")?;
    writeln!(&mut out_files[7], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino")?;
    writeln!(&mut out_files[8], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_qap, n_qap, m_qap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino")?;

    // Read graphs
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }
    let mut cthreads: Vec<JoinHandle<_>> = Vec::new();

    // Process graphs in different threads.
    for (graph, name) in graphs.clone() {
        cthreads.push(thread::spawn(move || {
            let mut dfvsi_org = DFVSInstance::new(graph.clone(), None, None);
            let start_time = Instant::now();
            dfvsi_org.compute_and_set_fast_upper(false);
            let time_upper = start_time.elapsed().as_millis();
            let start_time = Instant::now();
            dfvsi_org.compute_and_set_lower(false);
            let time_lower = start_time.elapsed().as_millis();
            (dfvsi_org, name, time_upper, time_lower)
    }));
    }


    let mut threads: Vec<Option<JoinHandle<Result<_,ThreadErr>>>> = Vec::new();
    let mut threads_info: Vec<(Receiver<_>, bool, usize, OsString)> = Vec::new();
    let mut timer: Vec<(Instant, Sender<u8>, bool)> = Vec::new();
    let mut instances = Vec::new();

    // Process graphs in different threads.
    for jh in cthreads {
        match jh.join() {
            Err(_) => {
                eprintln!("Some thread paniced");
            },
            Ok((dfvsi_org, name, upper_time, lower_time)) => {
                let upper_init = dfvsi_org.upper_bound.expect("was set");
                let lower_init = dfvsi_org.lower_bound.expect("was set");
                // Write general information 
                writeln!(out_files[0], "{:?}, {}, {}, {}, {}, {}, {}", name, 
                         dfvsi_org.graph.num_nodes(), 
                         dfvsi_org.graph.num_edges(), 
                         upper_init, upper_time, lower_init, lower_time)?;
                instances.push((dfvsi_org, name));
            },
        }
    }

    for (dfvsi_org, name) in instances {
        for p in 0..4 {
            let (start_sender, start_receiver) = channel();
            let (interrupt_sender, interrupt_receiver) = channel();
            let (done_sender, done_receiver) = channel();
            let n1 = name.clone();
            let mut dfvsi = dfvsi_org.clone();
            let priorities = priorities_org.clone();
            threads.push(Some(thread::spawn(move || {
                start_sender.send(1)?;
                if p == 0 || p == 1{
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[p], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            dfvsi.compute_and_set_fast_upper(true);
                            let upper = dfvsi.upper_bound.expect("was set");
                            eprintln!("Done {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((vec![dfvsi], vec![rule_stats], vec![upper])));
                        },
                        Err(_) => {
                            eprintln!("Interrupted {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(None);
                        },
                    };
                } else {
                    let mut kernels = Vec::new();
                    let mut rules = Vec::new();
                    let mut uppers = Vec::new();
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[2 + p-2], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            kernels.push(dfvsi.clone());
                            rules.push(rule_stats);
                        },
                        Err(_) => {
                            eprintln!("Interrupted {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((kernels, rules, uppers)));
                        },
                    };
                    // If `dfvsi` is reduced to zero stop here. 
                    if dfvsi.graph.num_nodes() == 0 {
                        uppers.push(dfvsi.solution.len());
                        eprintln!("Done {:?}, {}",n1, p);
                        done_sender.send(1)?;
                        return Ok(Some((kernels, rules, uppers)));
                    }
                    // Else compute heuristic ...
                    dfvsi.compute_and_set_fast_upper(true);
                    let upper = dfvsi.upper_bound.expect("was set");
                    uppers.push(upper);
                    // ... and continue with global lossy2 once + simple rules 
                    let global2 = dfvsi.apply_lossy_contract_globaly_once(2);
                    rules.push(vec![global2]);
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[0 + p-2], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            kernels.push(dfvsi.clone());
                            // resort rule_stats 
                            let mut new_stats: Vec<Option<RuleStats>> = vec![None;9];
                            for rs in rule_stats {
                                let pos;
                                match rs.rule {
                                    Rule::SimpleRules => pos = 0,
                                    Rule::Dome => pos = 1,
                                    Rule::SCC => pos = 2,
                                    Rule::AdvancedPetal => pos = 3,
                                    Rule::QuickAdvancedPetal => pos = 3,
                                    Rule::LinkNode => pos = 4,
                                    Rule::TwinNodes => pos = 5,
                                    Rule::Clique => pos = 6,
                                    Rule::Core => pos = 7,
                                    Rule::Dominion => pos = 8,
                                    _ => panic!("This rule should not be in here."),
                                }
                                new_stats[pos] = Some(rs.clone());
                            }
                            rules.push(new_stats.into_iter().map(|r| r.expect("all should have been filled")).collect());
                        },
                        Err(_) => {
                            eprintln!("Interrupted {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((kernels, rules, uppers)));
                        },
                    };
                    // If `dfvsi` is reduced to zero stop here. 
                    if dfvsi.graph.num_nodes() == 0 {
                        uppers.push(dfvsi.solution.len());
                        eprintln!("Done {:?}, {}",n1, p);
                        done_sender.send(1)?;
                        return Ok(Some((kernels, rules, uppers)));
                    }
                    // Else compute heuristic ...
                    dfvsi.compute_and_set_fast_upper(true);
                    let upper = dfvsi.upper_bound.expect("was set");
                    uppers.push(upper);
                    // ... and continue with lossy1 + simple rules 
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[2 + p-2], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            kernels.push(dfvsi.clone());
                            rules.push(rule_stats);
                            dfvsi.compute_and_set_fast_upper(true);
                            let upper = dfvsi.upper_bound.expect("was set");
                            uppers.push(upper);
                            eprintln!("Done {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((kernels, rules, uppers)));
                        },
                        Err(_) => {
                            eprintln!("Interrupted {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((kernels, rules, uppers)));
                        },
                    };
                }
            })));
            threads_info.push((done_receiver, false, p, name.clone()));

            // Try to join other threads until the current thread can start.
            let mut recvd = false;
            'outer: loop {
                // Check if message was received.
                if start_receiver.try_recv().is_ok() {
                    // put current time in vector 
                    recvd = true;
                    break 'outer
                }
                for i in 0..timer.len() {
                    let (recv, joined, go, g_name) = &mut threads_info[i];
                    if *joined {
                        continue;
                    }
                    // Join thread if it is finished:
                    if recv.try_recv().is_ok() {
                        let join_handle = threads[i].take().expect("`joined` is false");
                        match join_handle.join() {
                            Ok(Ok(Some((left_overs, rules, heur)))) => {
                                // Write file regarding of the `go` or how much was done before
                                // the interrupt.
                                if go == &0 || go == &1 {
                                    write_simple_stuff(dest, g_name, &left_overs[0], heur[0], &mut out_files[*go+1], &rules[0], *go)?;
                                } else {
                                    write_complex_stuff(dest, g_name, &left_overs, &heur, &mut out_files, &rules, *go)?;
                                }
                            },
                            Ok(Ok(None)) => {
                                if go == &0 || go == &1 {
                                    write_simple_empty(g_name, &mut out_files[*go+1])?;
                                } else {
                                    write_complex_stuff(dest, g_name, &vec![], &vec![], &mut out_files, &vec![], *go)?;
                                }
                            },
                            Ok(Err(_)) => eprintln!("Some thread paniced"),
                            Err(_) => eprintln!("Some thread paniced"),
                        }
                        *joined = true;
                        break 'outer
                    }
                    if !timer[i].2 && timer[i].0.elapsed() >= ultimate {
                        eprintln!("try to interrupt {}", i);
                        timer[i].1.send(1)?; 
                        timer[i].2 = true;
                    }
                }
            }
            if !recvd {
                eprintln!("trn to recv");
                start_receiver.recv()?;
            }
            // `start_receiver.try_recv()` should be ok by now.
            timer.push((Instant::now(), interrupt_sender, false));
        }
    }

    // All threads started, joining left overs:
    let mut remain = 1;
    while remain > 0 {
        remain = 0;
        for i in 0..timer.len() {
            let (recv, joined, go, name) = &mut threads_info[i];
            if *joined {
                continue;
            }
            if recv.try_recv().is_ok() {
                eprintln!("try to join {}", i);
                let join_handle = threads[i].take().expect("`joined` is false");
                match join_handle.join() {
                    Ok(Ok(Some((left_overs, rules, heur)))) => {
                        if go == &0 || go == &1 {
                            write_simple_stuff(dest, name, &left_overs[0], heur[0], &mut out_files[*go+1], &rules[0], *go)?;
                        } else {
                            write_complex_stuff(dest, name, &left_overs, &heur, &mut out_files, &rules, *go)?;
                        }
                    },
                    Ok(Ok(None)) => {
                        if go == &0 || go == &1 {
                            write_simple_empty(name, &mut out_files[*go+1])?;
                        } else {
                            write_complex_stuff(dest, name, &vec![], &vec![], &mut out_files, &vec![], *go)?;
                        }
                    },
                    Ok(Err(_)) => eprintln!("Some thread paniced"),
                    Err(_) => eprintln!("Some thread paniced"),
                }
                eprintln!("joined {}", i);
                *joined = true;
            } else if !timer[i].2 && timer[i].0.elapsed() >= ultimate {
                timer[i].1.send(1)?; 
                timer[i].2 = true;
            } else {
                remain += 1;
            }
        }
    }
    Ok(())
}

fn write_simple_stuff(
    dest: &str, g_name: &OsString, left_over: &DFVSInstance, heur: usize, out_file: &mut File, 
    rule_set: &Vec<RuleStats>, go: usize) -> Result<(), Box<dyn error::Error>> {
    left_over.graph.write_graph(File::create(format!("{}/{:?}_k_{}",dest,g_name, go))?)?;
    let mut line = String::new();
    line.push_str(&format!("{:?}, {}, {}, {}, {}, ", g_name, left_over.graph.num_nodes(), left_over.graph.num_edges(), left_over.solution.len(), heur));
    for r in 0..rule_set.len() {
        let rule = &rule_set[r];
        if r < rule_set.len()-1 {
            line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
        } else {
            line.push_str(&format!("{}, {}, {} ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
        }
    }
    writeln!(out_file, "{}",line)?;
    Ok(())
}

fn write_simple_empty(
    g_name: &OsString, out_file: &mut File) -> Result<(), Box<dyn error::Error>> {
    let mut line = String::new();
    line.push_str(&format!("{:?},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,", g_name));
    writeln!(out_file, "{}",line)?;
    Ok(())
}

fn write_complex_stuff(
    dest: &str, g_name: &OsString, left_overs: &Vec<DFVSInstance>, heurs: &Vec<usize>, 
    out_files: &mut Vec<File>, rule_set: &Vec<Vec<RuleStats>>, go: usize) -> Result<(), Box<dyn error::Error>> {
    let mut old_set = None;
    let mut rule_iter = rule_set.iter();
    for i in 0..3 {
        if left_overs.len() > i {
            left_overs[i].graph.write_graph(File::create(format!("{}/{:?}_k_{}",dest,g_name,go + i * 2))?)?;
            let mut line = String::new();
            line.push_str(&format!("{:?}, {}, {}, {}, {}, ", g_name, left_overs[i].graph.num_nodes(), left_overs[i].graph.num_edges(), left_overs[i].solution.len(), heurs[i]));
            // merge rule sets 
            if old_set.is_none() {
                old_set = Some(rule_iter.next().expect("`left_overs` still has elements").clone());
            } else {
                if i == 1 {
                    old_set = Some(RuleStats::merge_vecs(&old_set.expect("is some"), rule_iter.next().expect("`left_overs` still has elements"), true));
                }
                old_set = Some(RuleStats::merge_vecs(&old_set.expect("is some"), rule_iter.next().expect("`left_overs` still has elements"), true));
            }
            let merged_set = old_set.as_ref().expect("expect").clone();
            for r in 0..merged_set.len() {
                let rule = &merged_set[r];
                if r < merged_set.len()-1 {
                    if rule.rule == Rule::LossyClique(1) || 
                        rule.rule == Rule::GlobalLossyContraction(2) {
                        line.push_str(&format!("{}, {}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                    } else {
                        line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                    }
                } else {
                    if rule.rule == Rule::LossyClique(1) || 
                        rule.rule == Rule::GlobalLossyContraction(2) {
                        line.push_str(&format!("{}, {}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                    } else {
                        line.push_str(&format!("{}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                    }
                }
            }
            writeln!(out_files[go + 1 + i * 2], "{}",line)?;
        } else {
            let mut line = String::new();
            line.push_str(&format!("{:?} ", g_name));
            // file line with aproprate amount of ,
            match i {
                0 => line.push_str(",,,,,,,,,,,,,,,,,,,,"), 
                1 | 2=> line.push_str(",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"), 
                _ => (),
            }
            writeln!(out_files[go + 1 + i * 2], "{}",line)?;
        }
    }
    Ok(())
}
