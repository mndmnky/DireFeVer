use crate::digraph::{Digraph};
use crate::param_directed_feedback_arc_set::ParamDirectedFeedbackArcSetInstance;
use fxhash::{FxHashSet};
use std::collections::HashSet;

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ParamDirectedFeedbackVertexSetInstance {
    pub (crate) graph: Digraph,
    // Intermediate solution, contains nodes that are not any longer contained in `graph`.
    pub (crate) inter_solution: FxHashSet<usize>,
    // Parameter `k` does not include the elements in `inter_solution`.
    pub (crate) k: usize,
}

impl ParamDirectedFeedbackVertexSetInstance {

    /// Returns a new instance of the parametrized directed feedback set problem.
    pub fn new(graph: Digraph, k: usize) -> Self {
        ParamDirectedFeedbackVertexSetInstance {
            graph,
            inter_solution: FxHashSet::default(),
            k,
        }
    }

    /// Solves the directed feedback vertex set problem by transforming to a directed feedback arc
    /// set instance and solving that by iterative compression.
    pub fn solve_by_arc_version(&self) -> Option<FxHashSet<usize>> {
        let node_threshold = self.graph.num_nodes();
        let transformed: ParamDirectedFeedbackArcSetInstance = self.clone().into();
        transformed.solve().map(|sol| ParamDirectedFeedbackVertexSetInstance::transform_arc_solution(node_threshold, &sol))
    }

    /// Transforms a directed feedback vertex set instance to a directed feedback arc set
    /// instance. 
    fn transform_arc_solution(node_threshold: usize, arc_solution: &HashSet<(usize, usize)>) 
        -> FxHashSet<usize> {
            arc_solution.iter().map(|(src, trg)| {
                if *trg >= node_threshold {
                    *src
                } else {
                    *trg
                }
            }).collect::<FxHashSet<usize>>()
        }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::digraph::Digraph;
    use crate::param_directed_feedback_arc_set::ParamDirectedFeedbackArcSetInstance;

    #[test]
    fn transform_test() {
        let gr = Cursor::new("9 15 0\n3\n1 4\n2 5\n3\n4 6\n5 7 9\n9\n6\n8 6\n");
        let gr_trg = Cursor::new("18 24 0\n10\n11\n12\n13\n14\n15\n16\n\
                                17\n18\n3\n1 4\n2 5\n3\n4 6\n\
                                5 7 9\n9\n6\n8 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let g_trg = Digraph::read_graph(gr_trg);
        assert!(g_trg.is_ok());
        let g_trg = g_trg.unwrap();
        let dfvsi = ParamDirectedFeedbackVertexSetInstance::new(g, 2);
        let dfasi: ParamDirectedFeedbackArcSetInstance = dfvsi.into();
        assert_eq!(dfasi.graph, g_trg);
    }

    #[test]
    fn param_directed_feedback_vertex_set_test() {
        let gr = Cursor::new("9 15 0\n3\n1 4\n2 5\n3\n4 6\n5 7 9\n9\n6\n8 6\n");
        let g = Digraph::read_graph(gr);
        assert!(g.is_ok());
        let g = g.unwrap();
        let instance = ParamDirectedFeedbackVertexSetInstance::new(g, 2);
        let sol = instance.solve_by_arc_version();
        assert_eq!(sol, Some(vec![2usize, 5].into_iter().collect::<FxHashSet<usize>>()));
    }
}
