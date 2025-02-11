use crate::concepts_protocols::arithmetic_gate::gate::{Gate, Operation};
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
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

    pub fn evaluate(&self, inputs: Vec<T>) -> Vec<MultiLinearPolynomial<T>> {
        let mut evaluation_layers = vec![MultiLinearPolynomial::new(inputs.clone())];
        let mut running_inputs = inputs;

        self.layers.iter().for_each(|gates| {
            let mut next_inputs: Vec<T> = vec![];

            gates.iter().for_each(|gate| {
                let output = match gate.operation {
                    Operation::Add => running_inputs[gate.left] + running_inputs[gate.right],
                    Operation::Mul => running_inputs[gate.left] * running_inputs[gate.right],
                };

                next_inputs.push(output);
            });

            evaluation_layers.push(MultiLinearPolynomial::new(next_inputs.clone()));
            running_inputs = next_inputs;
        });

        println!("Evaluation of layers {:?}", evaluation_layers);
        evaluation_layers
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

        assert_eq!(
            *eval_output.last().unwrap().get_evaluation_points(),
            vec![Fq::from(15)]
        );
    }
}
