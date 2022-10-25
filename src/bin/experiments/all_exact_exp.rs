//!
//! This binary is only meant to run experiments and record their results.
//! For each input graph, several output files are created:
//! 1. One for each exact rule together with the simple rule.
//! 2. One with and one without the core rule
//! 3. One with and one without the advanced petal rule
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
        vec![Rule::SimpleRules],
        vec![Rule::SimpleRules, Rule::LinkNode],
        vec![Rule::SimpleRules, Rule::TwinNodes],
        vec![Rule::SimpleRules, Rule::Dome],
        vec![Rule::SimpleRules, Rule::Clique],
        vec![Rule::SimpleRules, Rule::Core],
        vec![Rule::SimpleRules, Rule::Dominion],
        vec![Rule::SimpleRules, Rule::SCC],
        vec![Rule::SimpleRules, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::QuickAdvancedPetal],
        vec![Rule::SimpleRules, Rule::Petal],
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Dominion, Rule::SCC, Rule::AdvancedPetal],
        vec![Rule::SimpleRules, Rule::LinkNode, Rule::TwinNodes, Rule::Dome, Rule::Clique, Rule::Core, Rule::Dominion, Rule::SCC],
    ];

    // Initialize output files
    let mut out_files = vec![
        File::create(format!("{}/sr.csv",dest))?,
        File::create(format!("{}/link.csv",dest))?,
        File::create(format!("{}/twins.csv",dest))?,
        File::create(format!("{}/dome.csv",dest))?,
        File::create(format!("{}/clique.csv",dest))?,
        File::create(format!("{}/core.csv",dest))?,
        File::create(format!("{}/dominion.csv",dest))?,
        File::create(format!("{}/scc.csv",dest))?,
        File::create(format!("{}/ap.csv",dest))?,
        File::create(format!("{}/qap.csv",dest))?,
        File::create(format!("{}/petal.csv",dest))?,
        File::create(format!("{}/all_but_core.csv",dest))?,
        File::create(format!("{}/all_but_petal.csv",dest))?,
    ];
    writeln!(&mut out_files[0], "name, nk, mk, sk,\
             t_st, n_st, m_st")?;

    writeln!(&mut out_files[1], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln")?;

    writeln!(&mut out_files[2], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_tn, n_tn, m_tn")?;

    writeln!(&mut out_files[3], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_dome, n_dome, m_dome")?;

    writeln!(&mut out_files[4], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_cliq, n_cliq, m_cliq")?;

    writeln!(&mut out_files[5], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_core, n_core, m_core")?;

    writeln!(&mut out_files[6], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_domino, n_domino, m_domino")?;

    writeln!(&mut out_files[7], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_scc, n_scc, m_scc")?;

    writeln!(&mut out_files[8], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_ap, n_ap, m_ap")?;
    
    writeln!(&mut out_files[9], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_qap, n_qap, m_qap")?;

    writeln!(&mut out_files[10], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_p, n_p, m_p")?;

    writeln!(&mut out_files[11], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc,\
             t_ap, n_ap, m_ap")?;

    writeln!(&mut out_files[12], "name, nk, mk, sk,\
             t_st, n_st, m_st,\
             t_ln, n_ln, m_ln,\
             t_tn, n_tn, m_tn,\
             t_dome, n_dome, m_dome,\
             t_cliq, n_cliq, m_cliq,\
             t_core, n_core, m_core,\
             t_domino, n_domino, m_domino,\
             t_scc, n_scc, m_scc")?;

    // Read graphs
    let mut graphs = Vec::new();
    for file in files {
        let graph = Digraph::read_graph(BufReader::new(File::open(file.clone())?))?;
        let name = file.file_stem().expect("Not a file.");
        graphs.push((graph, name.to_owned()));
    }
    let mut cthreads: Vec<JoinHandle<_>> = Vec::new();

    // Bounds not needed here
    // Process graphs in different threads.
    for (graph, name) in graphs.clone() {
        cthreads.push(thread::spawn(move || {
            let dfvsi_org = DFVSInstance::new(graph.clone(), None, None);
            (dfvsi_org, name)
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
            Ok((dfvsi_org, name)) => {
                instances.push((dfvsi_org, name));
            },
        }
    }

    for (dfvsi_org, name) in instances {
        for p in 0..13 {
            let (start_sender, start_receiver) = channel();
            let (interrupt_sender, interrupt_receiver) = channel();
            let (done_sender, done_receiver) = channel();
            let n1 = name.clone();
            let mut dfvsi = dfvsi_org.clone();
            let priorities = priorities_org.clone();
            threads.push(Some(thread::spawn(move || {
                start_sender.send(1)?;
                match dfvsi.exhaustive_fine_rules_stats(&priorities[p], &interrupt_receiver) {
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
                            Ok(Ok(Some((left_overs, rules)))) => {
                                // Write file regarding of the `go` or how much was done before
                                // the interrupt.
                                write_simple_stuff(dest, g_name, &left_overs, &mut out_files[*go], &rules, *go)?;
                            },
                            Ok(Ok(None)) => {
                                write_simple_empty(g_name, &mut out_files[*go], *go)?;
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
                    Ok(Ok(Some((left_overs, rules)))) => {
                        // Write file regarding of the `go` or how much was done before
                        // the interrupt.
                        write_simple_stuff(dest, name, &left_overs, &mut out_files[*go], &rules, *go)?;
                    },
                    Ok(Ok(None)) => {
                        write_simple_empty(name, &mut out_files[*go], *go)?;
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

fn write_simple_stuff(
    dest: &str, g_name: &OsString, left_over: &DFVSInstance, out_file: &mut File, 
    rule_set: &Vec<RuleStats>, go: usize) -> Result<(), Box<dyn error::Error>> {
    left_over.graph.write_graph(File::create(format!("{}/{:?}_k_{}",dest,g_name, go))?)?;
    let mut line = String::new();
    line.push_str(&format!("{:?}, {}, {}, {}, ", g_name, left_over.graph.num_nodes(), left_over.graph.num_edges(), left_over.solution.len()));
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
    g_name: &OsString, out_file: &mut File, go: usize) -> Result<(), Box<dyn error::Error>> {
    let mut line = String::new();
    if go < 10 {
        if go == 0 {
            line.push_str(&format!("{:?},,,,,,", g_name));
        } else {
            line.push_str(&format!("{:?},,,,,,,,,", g_name));
        }
    } else {
        line.push_str(&format!("{:?},,,,,,,,,,,,,,,,,,,,,,,,,,,", g_name));
    }
    writeln!(out_file, "{}",line)?;
    Ok(())
}

