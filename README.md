# DireFeVer

## a DIREcted FEedback VERtex set solver

DireFeVer is a solver for the directed vertex feedback set problem. This library provides the following algorithms and data structures:

### Data structures
* A directed, simple graph data structure that is implemented with two adjacency sets for each node.
* A directed feedback vertex set instance, that, beyond others, allows for redoing of alterations on the underlying graph data structure.

### Kernalization
* Simple kernalization rules as partially described in [^fn2] that include:
	* Removing multi-edges.
	* Adding loop nodes to the solution and removing them from the graph.
	* Removing sink and source nodes.
	* Merging nodes with incomming, or outgoing degree of 1 into its sole in-(resp. out-)neighbor. 

* Strongly connected component rules, that include:
	* Removal of edges connecting strongly connected components.
	* Removal of edges that connect strongly connected components in a subgraph, where strong edges were removed. This rule is described in [^fn3].

* `k`-flower rules, that include:
	* Variations of the removal of nodes that are connected to `k`+1 node disjunct cycles. Where `k` denotes a upper bound of the DFVS instance. This rule is described in [^fn4].
	* Variations of the removal of nodes that are connected to `k`+1-`l` node disjunct cycles `C`. Where `k` denotes a upper bound of the DFVS instance, and `l` denotes a lower bound of the subgraph, after the deletion of all nodes in `C`.

* Core rules, that include:
	* Variations of the clique rule, as described in [^fn3], there called "core" rule.
	* Variations of our core rule, which is a generalized version of the clique rule.

* The dome rule, as described in [^fn3].

### Heuristics
* Lower bound heuristics:
	* A heuristic that counts and removes small cycles.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| - 1 to the lower bound. This heuristic is described in [^fn3].

* Upper bound heuristics:
	* Variations of the big degree heuristic.
	* The lower bound heuristics that finds a clique `C`, removes it, and adds |`C`| to the upper bound.

### Exact Branching Algorithms
* A simple branching algorithm that should not currently be used.
* An advanced branching algorithm that first branches on cliques, then on daisies, and finally on small cycles.

### Parameterized Algorithm
* A solver that transforms the instance into a directed arc feedback set instance, then solves by iterative compression and transformation into skew edge multicut instances. As described by [^fn1]

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
* Further speed up of the scc- and advanced scc rule.

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
* binary for exact statistics

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
* SCC split + vc heursitic in exact bin

### 1.4.1
* Via vertex_cover now also computes an upper bound

### 1.5.0 (first open release)
* Removed `RebuildGraph`.
* All rebuild options are now in `DFVSInstance`.
* Removed all statistics.
* Removed the possibility to interrupt.

## Todo

### Some time
* Interrupter in other functions.
* Think about what priority rule list to use when.
* bound search space of bst if current_sol + lower_bound > best solution -> return best solution
* Better compute and set best bounds.
* Polished branch and reduce algorithm
* Polish and clean added stuff
* Include all remaining heuristics to the heuristic binary.
* Better structure for DFVS_instance and branch_instance.
* Rules from poly kernel with udfvs solution for heuristic solution.
* Test for heuristics
* Single instance heuristic stats need time and lower bound
* Interrupt does not seem to work for BST


[^fn1]: Chen, Jianer, et al. "A fixed-parameter algorithm for the directed feedback vertex set problem." Proceedings of the fortieth annual ACM symposium on Theory of computing. 2008.

[^fn2]: Levy, Hanoch, and David W. Low. "A contraction algorithm for finding small cycle cutsets." Journal of algorithms 9.4 (1988): 470-493.

[^fn3]: Lin, Hen-Ming, and Jing-Yang Jou. "On computing the minimum feedback vertex set of a directed graph by contraction operations." IEEE Transactions on computer-aided design of integrated circuits and systems 19.3 (2000): 295-307.

[^fn4]: Fleischer, Rudolf, Xi Wu, and Liwei Yuan. "Experimental study of FPT algorithms for the directed feedback vertex set problem." European Symposium on Algorithms. Springer, Berlin, Heidelberg, 2009.
