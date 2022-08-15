//!
//! This binary is only meant for experiments.
//! For each input graph, several output files are created:
//! 1. A csv containing statistics of the essential, the lossy clique and the lossy cycle rule and the resulting kernel for each input graph.
//! 2. A csv containing statistics of all the rules in step 1. and a single application of the
//!    lossy cut rule with complete exact rules afterwards, and the resulting kernel for each
//!    kernel from 1. that could not have been completely solved.
//! 3. Same as 2. but using the lossy indie cycle instead of the lossy cut rule
//! 4. Same as 2. but using the lossy semi indie cycle instead of the lossy cut rule
//! 5. With the basis of 2., 3. or 4. we apply the lossy contraction rule once with all exact rules
//!    afterwards.
//! 6. With the basis of 2., 3. or 4. we apply the lossy merge rule once with all exact rules
//!    afterwards.

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
    let m = Command::new("experiments_v2")
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
        vec![Rule::SimpleRules, Rule::LossyClique(1), Rule::LossyCycle(3) Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
    ];

    // TODO: from here

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/type_c_lossy.csv",dest))?,
        File::create(format!("{}/cut_lossy.csv",dest))?,
        File::create(format!("{}/indie_lossy.csv",dest))?,
        File::create(format!("{}/semi_indie_lossy.csv",dest))?,
        File::create(format!("{}/cut_lossy.csv",dest))?,
        File::create(format!("{}/indie_lossy_ffdsf.csv",dest))?,
        File::create(format!("{}/semi_indie_lossy.csv",dest))?,
        File::create(format!("{}/cut_lossy.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;
    writeln!(&mut out_files[1], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_ap, n_ap, m_ap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2")?;
    writeln!(&mut out_files[2], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lossy1, n_lossy1, m_lossy1, maxoff_lossy1,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_ap, n_ap, m_ap,\
             t_1lossy2, n_1lossy2, m_1lossy2, maxoff_1lossy2")?;

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
            (dfvsi_org, name)
    }));
    }


    let mut threads: Vec<Option<JoinHandle<Result<_,ThreadErr>>>> = Vec::new();
    let mut threads_info: Vec<(Receiver<_>, bool, OsString)> = Vec::new();
    let mut timer: Vec<(Instant, Sender<u8>, bool)> = Vec::new();
    let mut instances = Vec::new();

    // Process graphs in different threads.
    for jh in cthreads {
        match jh.join() {
            Err(_) => {
                eprintln!("Some thread paniced");
            },
            Ok((dfvsi_org, name)) => {
                instances.push((dfvsi_org, name));
            },
        }
    }

    for (dfvsi_org, name) in instances {
        let (start_sender, start_receiver) = channel();
        let (interrupt_sender, interrupt_receiver) = channel();
        let (done_sender, done_receiver) = channel();
        let n1 = name.clone();
        let mut dfvsi = dfvsi_org.clone();
        let priorities = priorities_org.clone();
        threads.push(Some(thread::spawn(move || {
            start_sender.send(1)?;
            let mut kernels = Vec::new();
            let mut rules = Vec::new();
            let mut uppers = Vec::new();
            match dfvsi.exhaustive_fine_rules_stats(&priorities[1], &interrupt_receiver) {
                Ok(rule_stats) => {
                    kernels.push(dfvsi.clone());
                    rules.push(rule_stats);
                },
                Err(_) => {
                    eprintln!("Interrupted {:?}",n1);
                    done_sender.send(1)?;
                    return Ok(None);
                },
            };
            // Collect nodes in the solution.
            // let bonus = dfvsi.solution.len(); // Not currently used.
            // If `dfvsi` is reduced to zero stop here. 
            if dfvsi.graph.num_nodes() == 0 {
                uppers.push(dfvsi.solution.len());
                eprintln!("Done {:?}",n1);
                done_sender.send(1)?;
                return Ok(Some((kernels, rules, uppers)));
            }
            // Else compute heuristic ...
            dfvsi.compute_and_set_fast_upper(true);
            let upper = dfvsi.upper_bound.expect("was set");
            uppers.push(upper);
            // ... and continue with global lossy2 once + all rules 
            let global2 = dfvsi.apply_global_lossy2_once(2);
            rules.push(vec![global2]);
            match dfvsi.exhaustive_fine_rules_stats(&priorities[0], &interrupt_receiver) {
                Ok(rule_stats) => {
                    kernels.push(dfvsi.clone());
                    rules.push(rule_stats);
                },
                Err(_) => {
                    eprintln!("Interrupted {:?}",n1);
                    done_sender.send(1)?;
                    return Ok(Some((kernels, rules, uppers)));
                },
            };
            // Collect nodes in the solution.
            // let bonus = dfvsi.solution.len(); // Not currently used.
            // If `dfvsi` is reduced to zero stop here. 
            if dfvsi.graph.num_nodes() == 0 {
                uppers.push(dfvsi.solution.len());
                eprintln!("Done {:?}",n1);
                done_sender.send(1)?;
                return Ok(Some((kernels, rules, uppers)));
            }
            // Else compute heuristic ...
            dfvsi.compute_and_set_fast_upper(true);
            let upper = dfvsi.upper_bound.expect("was set");
            uppers.push(upper);
            // ... and continue with global lossy2 once + all rules 
            match dfvsi.exhaustive_fine_rules_stats(&priorities[1], &interrupt_receiver) {
                Ok(rule_stats) => {
                    kernels.push(dfvsi.clone());
                    rules.push(rule_stats);
                    dfvsi.compute_and_set_fast_upper(true);
                    let upper = dfvsi.upper_bound.expect("was set");
                    uppers.push(upper);
                    eprintln!("Done {:?}",n1);
                    done_sender.send(1)?;
                    return Ok(Some((kernels, rules, uppers)));
                },
                Err(_) => {
                    eprintln!("Interrupted {:?}",n1);
                    done_sender.send(1)?;
                    return Ok(Some((kernels, rules, uppers)));
                },
            };
        })));
        threads_info.push((done_receiver, false, name.clone()));

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
                let (recv, joined, g_name) = &mut threads_info[i];
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
                            write_complex_stuff(dest, g_name, &left_overs, &heur, &mut out_files, &rules)?;
                        },
                        Ok(Ok(None)) => {
                            write_complex_stuff(dest, g_name, &vec![], &vec![], &mut out_files, &vec![])?;
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

    // All threads started, joining left overs:
    let mut remain = 1;
    while remain > 0 {
        remain = 0;
        for i in 0..timer.len() {
            let (recv, joined, name) = &mut threads_info[i];
            if *joined {
                continue;
            }
            if recv.try_recv().is_ok() {
                eprintln!("try to join {}", i);
                let join_handle = threads[i].take().expect("`joined` is false");
                match join_handle.join() {
                    Ok(Ok(Some((left_overs, rules, heur)))) => {
                        write_complex_stuff(dest, name, &left_overs, &heur, &mut out_files, &rules)?;
                    },
                    Ok(Ok(None)) => {
                        write_complex_stuff(dest, name, &vec![], &vec![], &mut out_files, &vec![])?;
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

fn write_complex_stuff(
    dest: &str, g_name: &OsString, left_overs: &Vec<DFVSInstance>, heurs: &Vec<usize>, 
    out_files: &mut Vec<File>, rule_set: &Vec<Vec<RuleStats>>) -> Result<(), Box<dyn error::Error>> {
    let mut old_set = None;
    let mut rule_iter = rule_set.iter();
    for i in 0..3 {
        if left_overs.len() > i {
            left_overs[i].graph.write_graph(File::create(format!("{}/{:?}_k_{}",dest,g_name,i))?)?;
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
            writeln!(out_files[i], "{}",line)?;
        } else {
            let mut line = String::new();
            line.push_str(&format!("{:?}", g_name));
            // file line with aproprate amount of ,
            match i {
                // 5*4
                0 => line.push_str(",,,,,,,,,,,,,,,,,,,,"), 
                // 6 + 11*3
                1 | 2 => line.push_str(",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"), 
                _ => (),
            }
            writeln!(out_files[i], "{}",line)?;
        }
    }
    Ok(())
}
