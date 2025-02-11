use crate::concepts_protocols::arithmetic_gate::gate::{Gate, Operation};
use ark_ff::PrimeField;
use std::marker::PhantomData;

pub struct Circuit<T: PrimeField> {
    _marker: PhantomData<T>,
    layers: Vec<Vec<Gate>>,
}

impl<T: PrimeField> Circuit<T> {
    pub fn new(layers: Vec<Vec<Gate>>) -> Self {
        Self {
            _marker: PhantomData,
            layers,
        }
    }

    pub fn evaluate(&self, inputs: Vec<T>) -> Vec<T> {
        self.layers.iter().fold(inputs, |inputs, gates| {
            let mut next_inputs: Vec<T> = Vec::with_capacity(inputs.len());

            gates.iter().for_each(|gate| {
                let output = match gate.operation {
                    Operation::Add => inputs[gate.left] + inputs[gate.right],
                    Operation::Mul => inputs[gate.left] * inputs[gate.right],
                };

                next_inputs.push(output);
            });

            next_inputs
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn arithmetic_gate_test() {
        let circuit = Circuit::new(vec![
            vec![
                Gate::new(0, 1, Operation::Add),
                Gate::new(2, 3, Operation::Mul),
            ],
            vec![Gate::new(0, 1, Operation::Add)],
        ]);

        let eval_output =
            circuit.evaluate(vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)]);

        assert_eq!(eval_output, vec![Fq::from(15)]);
    }
}
