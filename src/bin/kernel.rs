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

use dfvs_solver::{digraph::Digraph,  dfvs_instance::DFVSInstance, reduction_rules::Rule};

pub fn main() -> Result<(), Box<dyn error::Error>> {
    let m = Command::new("statistics")
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
    let ultimate = Duration::from_secs(1800);
    let dest: &str = m.value_of("dest").unwrap();
    let priorities = vec![
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::Lossy(1), Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::Lossy(2), Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::AdvancedPetal, Rule::Lossy(1)]
    ];
    let mut out_files = vec![
        File::create(format!("{}/kern_rules.csv",dest))?,
        File::create(format!("{}/simp_rules_lossy1.csv",dest))?,
        File::create(format!("{}/simp_rules_lossy2.csv",dest))?,
        File::create(format!("{}/kern_rules_lossy1.csv",dest))?
    ];
    writeln!(&mut out_files[0], "name, n, m, nk, mk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[1], "name, n, m, nk, mk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[2], "name, n, m, nk, mk,\
             t_st, n_st, m_st,\
             t_lossy2, n_lossy2, m_lossy2, maxoff_lossy2,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[3], "name, n, m, nk, mk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1")?;
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }
    let mut threads: Vec<Option<JoinHandle<Result<_,SendError<u8>>>>> = Vec::new();
    let mut threads_info: Vec<(Receiver<_>, bool, usize, usize, usize, OsString)> = Vec::new();
    let mut timer: Vec<(Instant, Sender<u8>)> = Vec::new();

    for (graph, name) in graphs.clone() {
        for p in 0..priorities.len() {
            let gr = graph.clone();
            let priority = priorities[p].clone();
            let (start_sender, start_receiver) = channel();
            let (interrupt_sender, interrupt_receiver) = channel();
            let (done_sender, done_receiver) = channel();
            let n1 = name.clone();
            threads.push(Some(thread::spawn(move || {
                start_sender.send(1)?;
                eprintln!("Start {:?}, {}",n1, p);
                let mut dfvsi = DFVSInstance::new(gr.clone(), None, None);
                match dfvsi.exhaustive_rules_stats(&priority, interrupt_receiver) {
                    Ok(rule_stats) => {
                        eprintln!("Done {:?}, {}",n1, p);
                        done_sender.send(1)?;
                        return Ok(Some((dfvsi, rule_stats)));
                    },
                    Err(_) => {
                        eprintln!("Interrupted {:?}, {}",n1, p);
                        done_sender.send(1)?;
                        return Ok(None);
                    },
                };
            })));
            threads_info.push((done_receiver, false, p, graph.num_nodes(), graph.num_edges(), name.clone()));
            let mut recvd = false;
            'outer: loop {
                // Check if message was received.
                if start_receiver.try_recv().is_ok() {
                    // put current time in vector 
                    recvd = true;
                    break 'outer
                }
                for i in 0..timer.len() {
                    let (recv, joined, rule_set, n, m, g_name) = &mut threads_info[i];
                    if *joined {
                        continue;
                    }
                    if recv.try_recv().is_ok() {
                        let join_handle = threads[i].take().expect("`joined` is false");
                        match join_handle.join() {
                            Ok(Ok(Some((left_over, rule_stats)))) => {
                                let file = File::create(format!("{}/{:?}_k_{}.csv",dest,g_name,rule_set))?;
                                left_over.graph.write_graph(file)?;
                                let mut line = String::new();
                                line.push_str(&format!("{:?}, {}, {}, {}, {}, ",g_name, n, m, left_over.graph.num_nodes(), left_over.graph.num_edges()));
                                for r in 0..rule_stats.len() {
                                    let rule = &rule_stats[r];
                                    if r < rule_stats.len()-1 {
                                        if rule.rule == Rule::Lossy(1) || rule.rule == Rule::Lossy(2) {
                                            line.push_str(&format!("{}, {}, {}, {},",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                                        } else {
                                            line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                                        }
                                    } else {
                                        if rule.rule == Rule::Lossy(1) || rule.rule == Rule::Lossy(2) {
                                            line.push_str(&format!("{}, {}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                                        } else {
                                            line.push_str(&format!("{}, {}, {} ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                                        }
                                    }
                                }
                                writeln!(out_files[*rule_set], "{}",line)?;
                            },
                            Ok(Ok(None)) => {
                                let mut line = String::new();
                                line.push_str(&format!("{:?}, {}, {},,,,,,,,,,,,,,,,,,",g_name, n, m));
                                if rule_set == &0 || rule_set == &3 {
                                    line.push_str(",,,,,,,,,,,")
                                }
                                if rule_set == &3 {
                                    line.push_str(",,,,")
                                }
                                writeln!(out_files[*rule_set], "{}",line)?;
                            },
                            Ok(Err(_)) => eprintln!("Some thread paniced"),
                            Err(_) => eprintln!("Some thread paniced"),
                        }
                        *joined = true;
                        break 'outer
                    }
                    if timer[i].0.elapsed() >= ultimate {
                        timer[i].1.send(1)?; 
                        let join_handle = threads[i].take().expect("`joined` is false");
                        // interrupt join
                        join_handle.join().expect("should work").err();
                        let mut line = String::new();
                        line.push_str(&format!("{:?}, {}, {},,,,,,,,,,,,,,,,,,",g_name, n, m));
                        if rule_set == &0 || rule_set == &3 {
                            line.push_str(",,,,,,,,,,,")
                        }
                        if rule_set == &3 {
                            line.push_str(",,,,")
                        }
                        writeln!(out_files[*rule_set], "{}",line)?;
                        *joined = true;
                        break 'outer
                    }
                }
            }
            if !recvd {
                eprintln!("trn to recv");
                start_receiver.recv()?;
            }
            // `start_receiver.try_recv()` should be ok by now.
            timer.push((Instant::now(), interrupt_sender));
        }
    }

    let mut remain = 1;
    while remain > 0 {
        remain = 0;
        for i in 0..timer.len() {
            let (recv, joined, rule_set, n, m, name) = &mut threads_info[i];
            if *joined {
                continue;
            }
            if recv.try_recv().is_ok() {
                let join_handle = threads[i].take().expect("`joined` is false");
                match join_handle.join() {
                    Ok(Ok(Some((left_over, rule_stats)))) => {
                        left_over.graph.write_graph(File::create(format!("{}/{:?}_k_{}.csv",dest,name,rule_set))?)?;
                        let mut line = String::new();
                        line.push_str(&format!("{:?}, {}, {}, {}, {}, ",name, n, m, left_over.graph.num_nodes(), left_over.graph.num_edges()));
                        for r in 0..rule_stats.len() {
                            let rule = &rule_stats[r];
                            if r < rule_stats.len()-1 {
                                if rule.rule == Rule::Lossy(1) || rule.rule == Rule::Lossy(2) {
                                    line.push_str(&format!("{}, {}, {}, {},",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                                } else {
                                    line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                                }
                            } else {
                                if rule.rule == Rule::Lossy(1) || rule.rule == Rule::Lossy(2) {
                                    line.push_str(&format!("{}, {}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                                } else {
                                    line.push_str(&format!("{}, {}, {} ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                                }
                            }
                        }
                        writeln!(out_files[*rule_set], "{}",line)?;
                    },
                    Ok(Ok(None)) => {
                        let mut line = String::new();
                        line.push_str(&format!("{:?}, {}, {},,,,,,,,,,,,,,,,,,",name, n, m));
                        if rule_set == &0 || rule_set == &3 {
                            line.push_str(",,,,,,,,,,,")
                        }
                        if rule_set == &3 {
                            line.push_str(",,,,")
                        }
                        writeln!(out_files[*rule_set], "{}",line)?;
                    },
                    Ok(Err(_)) => eprintln!("Some thread paniced"),
                    Err(_) => eprintln!("Some thread paniced"),
                }
                *joined = true;
            } else if timer[i].0.elapsed() >= ultimate {
                timer[i].1.send(1)?; 
                let join_handle = threads[i].take().expect("`joined` is false");
                // is interrupted
                join_handle.join().expect("should work").err();
                let mut line = String::new();
                line.push_str(&format!("{:?}, {}, {},,,,,,,,,,,,,,,,,,",name, n, m));
                if rule_set == &0 || rule_set == &3 {
                    line.push_str(",,,,,,,,,,,")
                }
                if rule_set == &3 {
                    line.push_str(",,,,")
                }
                writeln!(out_files[*rule_set], "{}",line)?;
                *joined = true;
            } else {
                remain += 1;
            }
        }
    }
    Ok(())
}
