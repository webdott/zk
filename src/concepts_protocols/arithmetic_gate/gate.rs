#[derive(Debug)]
pub enum Operation {
    Add,
    Mul,
}

#[derive(Debug)]
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
