pub enum Operation {
    Add,
    Mul,
}

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
