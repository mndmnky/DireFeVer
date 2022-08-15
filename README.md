# DireFeVer

## a DIREcted FEedback VERtex set solver

DireFeVer is a solver for the directed vertex feedback set problem. This library provides the following algorithms and data structures:

### Data structures
* A directed, simple graph data structure that is implemented with two adjacency sets for each node.
* A directed feedback vertex set instance, that, beyond others, allows for redoing of alterations on the underlying graph data structure.

### Kernelization
* Simple kernelization rules as partially described in [^fn2] that include:
	* Removing multi-edges.
	* Adding loop nodes to the solution and removing them from the graph.
	* Removing sink and source nodes.
	* Merging nodes with incomming, or outgoing degree of 1 into its sole in-(resp. out-)neighbor. 

* Strongly connected component rules, that include:
	* Removal of edges connecting strongly connected components.
	* Removal of edges that connect strongly connected components in a subgraph, where strong edges were removed. This rule is described in [^fn3].

* Petal rules that, for a node `v`, counts the adjacent, otherwise node disjunct cycles, also called petals. The first two rules are described in [^fn4].
	* If `v` has one petal, `v` can be contracted.
	* If `v` has more petals than the effective upper bound, `v` can be added to the solution.
	* In an advanced version, if `v` has more petals than the effective upper bound minus the lower bound of the graph minus all nodes in `v`s petals (including) `v`, `v` can be added to the solution.

* Core rules, that include:
	* Variations of the clique rule, as described in [^fn3], there called "core" rule.
	* Variations of our core rule, which is a generalized version of the clique rule.

* The dome rule, as described in [^fn3].

* Vertex Cover rules applied on subgraphs with reciprocal edges:
	* The link node rule, that contract nodes adjacent to and only to exactly two reciprocal edges, going to and coming from nodes `v` and `w`. As described for the vertex cover problem in [^fn5]. To adapt this rule for the directed feedback vertex set problem, one can not use cases where there exists a directed path from `v` to `w` since that might lead to a directed circle that was not present before.
	* The twin node rule, as described in [^fn6]. This rule looks for a pair of nodes with identical 3 reciprocal neighbors `v`, `w` and `u`, if there exists a reciprocal edge between any two nodes of `v`, `w` and `u`, those three nodes can be added to the cover.
	* The unconfined (here dubbed dominion) rule as described in [^fn7].
	* And the crown rule as described in [^fn8].

* Lossy kernelization rules:
	* One, that adds all nodes within strong cliques of size 1+`quality`. The size of an optimal solution of the resulting kernel can become at worst 1 + 1/`quality` times as large.
	* A second one, which contracts nodes adjacent to at most `quality` petals. Each application of this rule can increase the solution by `quality`-1. If the rule is applied only once, on only cycle disjunct nodes, the size of an optimal solution of the resulting kernel can become at worst `quality` times as large. If other lossy rules follow this, however, the resulting solution can become `quality` times the approximation factor of the other rules as large as an optimal solution. This rule can be simplified by using nodes with a maximal incoming- or outgoing degree of `quality`. But this might miss some valid nodes.
	* TODO describe the other lossy rules!
	* These rules reset the current upper bound and current best solution. Or they should!

### Heuristics
* Lower bound heuristics:
	* A heuristic that counts and removes small cycles.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| - 1 to the lower bound. This heuristic is described in [^fn3].

* Upper bound heuristics:
	* Variations of the big degree heuristic.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| to the upper bound.

* Vertex-cover heuristic:
	* Transforms strongly connected components to a vertex cover instance by only regarding reciprocal edges as single undirected edges and discarding the rest. 
	* Then, we run a vertex cover solver over that instance. Here we use Duck and Cover [^fn9].
	* If the returned solution is also a solution for the original strongly connected component, then the solution is optimal.
	* Otherwise, we can use the solution as a lower bound, which can also be extended to an upper bound by solving the leftover graph.

### Exact Branching Algorithms
* An advanced branching algorithm that branches on the first available option in the following priority list: 
	1. A clique of size 3 or more.
	2. A link node on which the rule could not be applied.
	3. A daisy.
	4. The node with the highest chosen weight.

### Parameterized Algorithm
* A solver that transforms the instance into a directed arc feedback set instance, then solves by iterative compression and transformation into skew edge multicut instances. As described by [^fn1]

### Statistics
* Statistics for the reduction rules.

## Installation
* Make sure `cargo` is installed.
* Run `cargo build --release`.

## Changelog (at 1.0.0)

### 1.0.1 
* Speed up `apply_advanced_scc_rule()`.
* Added `apply_advanced_scc_rule()` to `BranchInstance`.
* Added `cai_weight()` and `lin_weight()` to `Digraph`.
* Added bottom-up and top-down weight heuristics.
* Added top-bottom-switch weight heuristics.
* Added local-search heuristics.
* Added clique upper bound branch heuristic.
* Simple statistics for the weight heuristics.
* Binary for statistics on the heuristics.

### 1.1.0
* Merge `DirectedFeedbackVertexSetInstance` and `BranchInstance`, only using `RebuildGraph`.
* Binary that records reduction statistics on multiple files.
* Iterative approach for scc- and advanced scc rule.
* Further, speed up of the scc- and advanced scc rule.

### 1.2.0 
* Added `exhaustive_reductions()`.
* Added `apply_exhaustive_core_rule()`.
* Added `compute_and_set_lower()`.
* Added `apply_local_k_daisy()`.
* Added `compute_and_set_simple_upper_lower()`.
* Added `exhaustive_local_search()`.
* Fixed `mult_reduction_stats.rs`
* Struct interrupter:
	* Only handles time_outs
	* Implemented for exhaustive reductions
* Overhauled `advanced_branching()`.
* Put `exhaustive_reductions()` in each heuristic and branch and bound. 
* Added `current_best` to `dfvs_instance`.
* Proper bin for statistics on single instance heuristics.
* `advanced_clique_heuristic()` branches only to a maximal depth of 3, which could still be too much.
* Fixed a bug where `{RebuildGraph,DFVSInstance}::_reset_reductions()` did not reset `changes` (resp. `reductions`).
* Changes to `advanced_branching()`: 
	* now branches on highest cai weight (0.3) node if no more daisy or clique remains.
	* checks if the current lower bound (+ the size of the current solution) is greater or equal than the current best solution. If so, returns the best current solution. 
	* Handles sccs separately.

### 1.2.1 (before interrupt in recursive branching)
* Added branching on double paths in `advanced_branching()`.

### 1.2.2
* Complete interrupt for BST.
* Binary for exact statistics.

### 1.3.0
* Link node rule + restructuring of DFVSInstance.

### 1.3.1
* Use fxhash for all hashmaps and hashsets that use usize as a key.

### 1.3.2
* Crown rule

### 1.3.3 
* Twin rule
* Dominion rule

### 1.4.0
* Changes when computing and updating the upper and lower bound.
* Heuristic over vertex cover
* SCC split + vc heuristic in exact bin

### 1.4.1
* Via vertex_cover now also computes an upper bound

### 1.5.0 (first open release)
* Removed `RebuildGraph`.
* All rebuild options are now in `DFVSInstance`.
* Removed all statistics.
* Removed the possibility to interrupt.

### 1.5.1 
* Updated README

### 1.6.0
* `count_patels()` implemented. Plus a version that returns the left-over graph.
* `apply_petal_rules()` and `apply_advanced_petal_rules()` implemented and integrated.
* `apply_lossy_rules()` implemented and integrated.
* `kernel` binary.
* `stats.rs`.

### 1.6.1 
* `kernel` now also records the current solution size.

### 1.7.0 
* `apply_simple_lossy2_rules()` and `apply_advanced_lossy2_rules()` implemented and integrated.
* Single node version of some reduction rules for better interruptability.
* Adapted `kernel` binary.

### 1.7.1 
* Added fast lossy to `kernel`

### 1.7.2 
* Added `experiments` bin with fixed experiments 
* Added `GlobalLossy2` rule

### 1.8.0
* Local simple rules
* Fast heuristic now actually fast. (Hopefully)

### 1.8.1
* Fast advanced petal rule.

### 1.8.2
* Different binaries for experiments.
* Some clean-up and tests.

### 1.9.0
* Find semi disjoint cycles 
* Lossy cycle rule (with supporting functions)
* Lossy indie cycle rule (with supporting functions)
* Lossy semi indie cycle rule (with supporting functions)
* Lossy merge rule (with supporting functions and DFVSInstance support)
* Binaries for experiments

### 1.9.1 
* Find indie cycles in BFS fashion instead.
* Added an advanced cut rule.
* And another cut variation.
* And cut with parameter.

## Todo

### 1.9.X
* Sub functions of the cut rule needs to be fixed (more efficient, remove redundant)
* working tests for both indie cycle rules
* Working test for adv cut rule 
* Quality check/guarentee for adv cut rule

### Very important but lots of work 
* Split sccs in while exhaustively applying the rules!

### Next
* Speed-up rules
* Twin-rule can merge nodes instead. Would take some work though.

### Some other things
* SCC in disjunct cycle heuristic can be simple. Check if simple is much faster than advanced.
* Twin nodes can also be merged if no edges are between their neighbors.
* Constant rule priority lists.
* UDFVS solution as heuristic.
* Poly kernel rules.

# Still needed?
* TEST `single versions`
* Check if BST can be bounded further (current solution + current best lower bound > best solution -> return best solution).
* Polish branch and reduce algorithm (can be contracted).

[^fn1]: Chen, Jianer, et al. "A fixed-parameter algorithm for the directed feedback vertex set problem." Proceedings of the fortieth annual ACM symposium on Theory of computing. 2008.

[^fn2]: Levy, Hanoch, and David W. Low. "A contraction algorithm for finding small cycle cutsets." Journal of algorithms 9.4 (1988): 470-493.

[^fn3]: Lin, Hen-Ming, and Jing-Yang Jou. "On computing the minimum feedback vertex set of a directed graph by contraction operations." IEEE Transactions on computer-aided design of integrated circuits and systems 19.3 (2000): 295-307.

[^fn4]: Fleischer, Rudolf, Xi Wu, and Liwei Yuan. "Experimental study of FPT algorithms for the directed feedback vertex set problem." European Symposium on Algorithms. Springer, Berlin, Heidelberg, 2009.

[^fn5]: Chen, Jianer, Iyad A. Kanj, and Weijia Jia. "Vertex cover: further observations and further improvements." Journal of Algorithms 41.2 (2001): 280-301.

[^fn6]: Xiao, Mingyu, and Hiroshi Nagamochi. "Confining sets and avoiding bottleneck cases: A simple maximum independent set algorithm in degree-3 graphs." Theoretical Computer Science 469 (2013): 92-104.

[^fn7]: Akiba, Takuya, and Yoichi Iwata. "Branch-and-reduce exponential/FPT algorithms in practice: A case study of vertex cover." Theoretical Computer Science 609 (2016): 211-225.

[^fn8]: Abu-Khzam, Faisal N., et al. "Kernelization algorithms for the vertex cover problem." (2017).

[^fn9]: https://github.com/mndmnky/duck-and-cover
