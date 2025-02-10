use ark_ff::PrimeField;

enum Operation {
    Add,
    Mul,
}

pub struct Gate {
    left: usize,
    right: usize,
    output: usize,
    operation: Operation,
}

pub struct Circuit<T: PrimeField> {
    layer_values: Vec<T>,
    layer_gates: Vec<Gate>,
}

impl<T: PrimeField> Circuit<T> {
    pub fn new(layer_values: Vec<T>, layer_gates: Vec<Gate>) -> Self {
        Self {
            layer_values,
            layer_gates,
        }
    }
}
