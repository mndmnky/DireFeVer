//!
//! This binary is purely meant for experiments.
//! For each input graph, several output files are created:
//! * A csv with general statistics of the graph. Such as number of nodes, number of edges, some
//! upper bound and some lower bound.
//! * A csv containing statistics of the essential rules and the resulting kernel.
//! * A csv containing statistics of the essential and the lossy1 rule and the resulting kernel.
//! * A csv containing statistics of the essential, the lossy1 rule and a single application of the
//! lossy2 rule without the lossy1 rule afterwards, and the resulting kernel.
//! * A csv containing statistics of the essential, the lossy1 rule and a single application of the
//! lossy2 rule, and the resulting kernel.
//! * A csv containing statistics of the essential, the lossy1 rule, a single application of the
//! lossy2 rule and a final heuristic on the remaining kernel which only removes as many nodes as
//! the essential rules added into the solution in the first place, and the resulting kernel.
//!
//! TODO: repeat with complete rules

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
        vec![Rule::SimpleRules, Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::Lossy(1), Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
    ];

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/init_stats.csv",dest))?,
        File::create(format!("{}/kern_rules.csv",dest))?,
        File::create(format!("{}/kern_rules_lossy1.csv",dest))?,
        File::create(format!("{}/kern_rules_lossy1_1lossy2.csv",dest))?,
        File::create(format!("{}/kern_rules_lossy1_1lossy2_4.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, n, m, upper_bound")?;
    writeln!(&mut out_files[1], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[2], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[3], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[4], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2")?;

    // Read graphs
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }

    let mut threads: Vec<Option<JoinHandle<Result<_,ThreadErr>>>> = Vec::new();
    let mut threads_info: Vec<(Receiver<_>, bool, usize, usize, usize, OsString)> = Vec::new();
    let mut timer: Vec<(Instant, Sender<u8>, bool)> = Vec::new();

    // Process graphs in different threads.
    for (graph, name) in graphs.clone() {
        let dfvsi_org = DFVSInstance::new(graph.clone(), None, None);
        dfvsi_org.compute_and_set_fast_upper(false);
        let upper_init = dfvsi_org.upper_bound.expect("was set");
        // Write general information 
        writeln!(out_files[0], "{:?}, {}, {}, {}", name, graph.num_nodes(), graph.num_edges(), upper_init)?;
        for p in 0..2 {
            let (start_sender, start_receiver) = channel();
            let (interrupt_sender, interrupt_receiver) = channel();
            let (done_sender, done_receiver) = channel();
            let n1 = name.clone();
            let mut dfvsi = dfvsi_org.clone();
            let priorities = priorities_org.clone();
            threads.push(Some(thread::spawn(move || {
                start_sender.send(1)?;
                if p == 0 {
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[0], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            let upper = dfvsi
                                .top_down_weight_heuristic(
                                    &Digraph::cai_weight, (0.2,0.0), &vec![Rule::SimpleRules], true
                                    );
                            eprintln!("Done {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(Some((vec![dfvsi], vec![rule_stats], vec![upper.len()])));
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
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[1], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            kernels.push(dfvsi.clone());
                            rules.push(rule_stats);
                        },
                        Err(_) => {
                            eprintln!("Interrupted {:?}, {}",n1, p);
                            done_sender.send(1)?;
                            return Ok(None);
                        },
                    };
                    // Collect nodes in the solution.
                    // let bonus = dfvsi.solution.len(); // Not currently used.
                    // If `dfvsi` is reduced to zero stop here. 
                    if dfvsi.graph.num_nodes() == 0 {
                        uppers.push(dfvsi.solution.len());
                        eprintln!("Done {:?}, {}",n1, p);
                        done_sender.send(1)?;
                        return Ok(Some((kernels, rules, uppers)));
                    }
                    // Else compute heuristic ...
                    let upper = dfvsi
                        .top_down_weight_heuristic(
                            &Digraph::cai_weight, (0.2,0.0), &vec![Rule::SimpleRules], true
                            );
                    uppers.push(upper.len());
                    // ... and continue with lossy1 + simple rules 
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[2], &interrupt_receiver) {
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
                    let upper = dfvsi
                        .top_down_weight_heuristic(
                            &Digraph::cai_weight, (0.2,0.0), &vec![Rule::SimpleRules], true
                            );
                    uppers.push(upper.len());
                    // ... and continue with global lossy2 once + simple rules 
                    let global2 = dfvsi.apply_global_lossy2_once(2);
                    rules.push(vec![global2]);
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[1], &interrupt_receiver) {
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
                    let upper = dfvsi
                        .top_down_weight_heuristic(
                            &Digraph::cai_weight, (0.2,0.0), &vec![Rule::SimpleRules], true
                            );
                    uppers.push(upper.len());
                    // ... and continue with lossy1 + simple rules 
                    match dfvsi.exhaustive_fine_rules_stats(&priorities[2], &interrupt_receiver) {
                        Ok(rule_stats) => {
                            kernels.push(dfvsi.clone());
                            rules.push(rule_stats);
                            let upper = dfvsi
                                .top_down_weight_heuristic(
                                    &Digraph::cai_weight, (0.2,0.0), &vec![Rule::SimpleRules], true
                                    );
                            uppers.push(upper.len());
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
            threads_info.push((done_receiver, false, p, graph.num_nodes(), graph.num_edges(), name.clone()));

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
                    let (recv, joined, go, _, _, g_name) = &mut threads_info[i];
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
                                if go == &0 {
                                    write_simple_stuff(dest, g_name, &left_overs[0], heur[0], &mut out_files[0], &rules[0])?;
                                } else {
                                    write_complex_stuff(dest, g_name, &left_overs, &heur, &mut out_files, &rules)?;
                                }
                            },
                            Ok(Ok(None)) => {
                                if go == &0 {
                                    write_simple_empty(g_name, &mut out_files[0])?;
                                } else {
                                    write_complex_stuff(dest, g_name, &vec![], &vec![], &mut out_files, &vec![])?;
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
            let (recv, joined, go, _, _, name) = &mut threads_info[i];
            if *joined {
                continue;
            }
            if recv.try_recv().is_ok() {
                eprintln!("try to join {}", i);
                let join_handle = threads[i].take().expect("`joined` is false");
                match join_handle.join() {
                    Ok(Ok(Some((left_overs, rules, heur)))) => {
                        if go == &0 {
                            write_simple_stuff(dest, name, &left_overs[0], heur[0], &mut out_files[0], &rules[0])?;
                        } else {
                            write_complex_stuff(dest, name, &left_overs, &heur, &mut out_files, &rules)?;
                        }
                    },
                    Ok(Ok(None)) => {
                        if go == &0 {
                            write_simple_empty(name, &mut out_files[0])?;
                        } else {
                            write_complex_stuff(dest, name, &vec![], &vec![], &mut out_files, &vec![])?;
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
    rule_set: &Vec<RuleStats>) -> Result<(), Box<dyn error::Error>> {
    left_over.graph.write_graph(File::create(format!("{}/{:?}_k_0",dest,g_name))?)?;
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
    // 9 * 3 - 1
    line.push_str(&format!("{:?},,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,", g_name));
    writeln!(out_file, "{}",line)?;
    Ok(())
}

fn write_complex_stuff(
    dest: &str, g_name: &OsString, left_overs: &Vec<DFVSInstance>, heurs: &Vec<usize>, 
    out_files: &mut Vec<File>, rule_set: &Vec<Vec<RuleStats>>) -> Result<(), Box<dyn error::Error>> {
    let mut old_set = None;
    for i in 0..4 {
        if left_overs.len() > i {
            left_overs[i].graph.write_graph(File::create(format!("{}/{:?}_k_{}",dest,g_name,i+1))?)?;
            let mut line = String::new();
            line.push_str(&format!("{:?}, {}, {}, {}, {}, ", g_name, left_overs[i].graph.num_nodes(), left_overs[i].graph.num_edges(), left_overs[i].solution.len(), heurs[i]));
            // merge rule sets 
            if old_set.is_none() {
                old_set = Some(rule_set[i].clone());
            } else {
                let order = i == 2; 
                old_set = Some(RuleStats::merge_vecs(&old_set.expect("is some"), &rule_set[i].clone(), order));
            }
            for r in 0..rule_set[i].len() {
                let rule = &rule_set[i][r];
                if r < rule_set[i].len()-1 {
                    if rule.rule == Rule::Lossy(1) || 
                        rule.rule == Rule::Lossy(2) || 
                            rule.rule == Rule::SimpleLossy2(2) || 
                            rule.rule == Rule::AdvancedLossy2(2) {
                        line.push_str(&format!("{}, {}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                    } else {
                        line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                    }
                } else {
                    if rule.rule == Rule::Lossy(1) || 
                        rule.rule == Rule::Lossy(2) || 
                            rule.rule == Rule::SimpleLossy2(2) || 
                            rule.rule == Rule::AdvancedLossy2(2) {
                        line.push_str(&format!("{}, {}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                    } else {
                        line.push_str(&format!("{}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                    }
                }
                writeln!(out_files[i], "{}",line)?;
            }
        } else {
            let mut line = String::new();
            line.push_str(&format!("{:?}, ", g_name));
            // file line with aproprate amount of ,
            match i {
                // 4 * 3 - 1
                0 => line.push_str(",,,,,,,,,,,,,,,"), //??
                // 4 * 3 - 1 + 4
                1 => line.push_str(",,,,,,,,,,,,,,,,,,,"), //??
                // 4 * 3 - 1 + 4*2
                2 | 3 => line.push_str(",,,,,,,,,,,,,,,,,,,,,,,"), //??
                _ => (),
            }
            writeln!(out_files[i], "{}",line)?;
        }
    }
    Ok(())
}
