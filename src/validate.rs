//! 
//! A module for simple validation of a dfvs solution S, by checking if the removal of each node in
//! S results into a DAG.

use crate::dfvs_instance::DFVSInstance;
use fxhash::FxHashSet;

impl DFVSInstance {

    /// Validates a solution `sol` for a untouched DFVSInstance. 
    /// Does not check for optimality. 
    ///
    /// # Panics 
    /// Panics if `self` already holds elements in `self.solution`.
    pub fn validate(&self, sol: &FxHashSet<usize>) -> bool {
        assert!(self.solution.is_empty());
        let mut c_ins = self.clone();
        c_ins.graph.remove_nodes(sol.iter().copied());
        c_ins.apply_simple_rules();
        let mut is_solution = c_ins.solution.is_empty(); 
        is_solution = is_solution && c_ins.graph.num_nodes() == 0;
        is_solution
    }

    /// Validates a solution `sol` for the left over `graph` of DFVSInstance. 
    /// Does not check for optimality. 
    pub fn validate_left_over(&self, sol: &FxHashSet<usize>) -> bool {
        let mut c_ins = DFVSInstance::new(self.graph.clone(), None, None);
        c_ins.graph.remove_nodes(sol.iter().copied());
        c_ins.apply_simple_rules();
        let mut is_solution = c_ins.solution.is_empty(); 
        is_solution = is_solution && c_ins.graph.num_nodes() == 0;
        eprintln!("is solution {:?}", is_solution);
        is_solution
    }

}
