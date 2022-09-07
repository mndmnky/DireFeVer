//!
//! This binary is only meant for experiments.
//! Runs Type C lossy reduction rules and outputs a csv containing statistics for the rules and the
//! resulting kernel.

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
    let m = Command::new("typec")
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

    // Rule priority sets
    let priorities_org = vec![
        //vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::LossyClique(1), Rule::LossyCycle(3), Rule::Dome, Rule::SCC, Rule::AdvancedPetal],
    ];

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/type_c_rules.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, nk, mk, sk, uk,\
             t_st, n_st, m_st,\
             t_lcli, n_lcli, m_lcli, maxoff_lcli,\
             t_lcy, n_lcy, m_lcy, maxoff_lcy,\
             t_dome, n_dome, m_dome,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;

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
            let dfvsi_org = DFVSInstance::new(graph.clone(), None, None);
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
            match dfvsi.exhaustive_fine_rules_stats(&priorities[0], &interrupt_receiver) {
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
                    return Ok(None);
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
                remain += 1;
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
    if rule_set.len() > 0{
        left_overs[0].graph.write_graph(File::create(format!("{}/{:?}_t1k",dest,g_name))?)?;
        let mut line = String::new();
        line.push_str(&format!("{:?}, {}, {}, {}, {}, ", g_name, left_overs[0].graph.num_nodes(), left_overs[0].graph.num_edges(), left_overs[0].solution.len(), heurs[0]));
        for r in 0..rule_set[0].len() {
            let rule = &rule_set[0][r];
            if r < rule_set[0].len()-1 {
                if rule.rule == Rule::LossyClique(1) || 
                    rule.rule == Rule::LossyCycle(3) {
                    line.push_str(&format!("{}, {}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                } else {
                    line.push_str(&format!("{}, {}, {}, ",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                }
            } else {
                if rule.rule == Rule::LossyClique(1) || 
                    rule.rule == Rule::LossyCycle(3) {
                    line.push_str(&format!("{}, {}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges, rule.suc_apps));
                } else {
                    line.push_str(&format!("{}, {}, {}",rule.time_took, rule.reduced_nodes, rule.reduced_edges));
                }
            }
        }
        writeln!(out_files[0], "{}",line)?;
    } else {
        let mut line = String::new();
        line.push_str(&format!("{:?}", g_name));
        // file line with aproprate amount of ,
        // 6*4
        line.push_str(",,,,,,,,,,,,,,,,,,,,,,,,");
        writeln!(out_files[0], "{}",line)?;
    }
    Ok(())
}
