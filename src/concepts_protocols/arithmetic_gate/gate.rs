#[derive(Debug, Clone)]
pub enum Operation {
    Add,
    Mul,
}

#[derive(Debug, Clone)]
pub struct Gate {
    pub left: usize,
    pub right: usize,
    pub operation: Operation,
}

impl Gate {
    pub fn new(left: usize, right: usize, operation: Operation) -> Self {
        Self {
            left,
            right,
            operation,
        }
    }
}
