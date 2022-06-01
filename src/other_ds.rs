use std::iter::FromIterator;

/// Stores a set of nodes.
///
/// Implemented with a `Vec<usize>`. Trade-offs apply.
#[derive(Debug, Default, Clone)]
pub struct NodeSet {
    set: Vec<bool>,
    len: usize,
}

impl NodeSet {

    /// Returns a new node set.
    pub fn new() -> Self {
        Self {
            set: Vec::new(),
            len: 0,
        }
    }

    /// Creats a nodes set with an initial capacity of `capacity`.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            set: Vec::with_capacity(capacity),
            len: 0,
        }
    }

    /// Check if `self` contains `node`.
    pub fn contains(&self, &node: &usize) -> bool {
        if node >= self.set.len() {
            false
        } else {
            self.set[node]
        }
    }

    /// Returns the number of elements in the set.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the set contains no elements.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Inserts `node` into the set.
    pub fn insert(&mut self, node: usize) {
        if node < self.set.len() {
            if !self.set[node] {
                self.set[node] = true;
                self.len += 1;
            }
        } else {
            self.set.resize(node+1, false);
            self.set[node] = true;
            self.len += 1;
        }
    }

	/// Removes `node` from the set.
	pub fn remove(&mut self, node: usize) {
        if node < self.set.len() && self.set[node] {
            self.set[node] = false;
            self.len -= 1;
        }
	}

}

impl IntoIterator for NodeSet {
    type Item = usize;
    type IntoIter = NodeSetIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        NodeSetIntoIter {
            set: self.set,
            len: self.len,
            i: 0,
        }
    }

}

impl FromIterator<usize> for NodeSet {

    fn from_iter<I: IntoIterator<Item=usize>>(iter: I) -> Self {
        let mut c = NodeSet {
            set: Vec::new(),
            len: 0,
        };
        for i in iter {
            c.insert(i);
        }
        c
    }

}

pub struct NodeSetIntoIter {
    set: Vec<bool>,
    len: usize,
    i: usize
}

impl Iterator for NodeSetIntoIter {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.len == 0 { return None; }
        self.len -= 1;
        while !self.set[self.i] {
            self.i += 1;
        }
        self.i += 1;
        Some(self.i - 1)
    }

}
