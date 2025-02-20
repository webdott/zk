use crate::concepts_protocols::arithmetic_gate::gate::{Gate, Operation};
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use ark_ff::PrimeField;
use std::cmp::max;

use std::marker::PhantomData;

#[derive(Clone)]
pub struct Circuit<T: PrimeField> {
    _marker: PhantomData<T>,
    layers: Vec<Vec<Gate>>,
    layer_evaluations: Vec<MultiLinearPolynomial<T>>,
}

impl<T: PrimeField> Circuit<T> {
    pub fn new(layers: Vec<Vec<Gate>>) -> Self {
        Self {
            _marker: PhantomData,
            layers,
            layer_evaluations: vec![],
        }
    }

    pub fn get_layer_count(&self) -> usize {
        self.layers.len()
    }

    pub fn evaluate_at_input(&mut self, inputs: Vec<T>) -> MultiLinearPolynomial<T> {
        let mut evaluation_layers = vec![MultiLinearPolynomial::new(inputs.clone())];
        let mut running_inputs = inputs;

        self.layers.iter().for_each(|gates| {
            let mut next_inputs: Vec<T> =
                vec![T::from(0); max(gates.len().next_power_of_two(), 2) as usize];

            gates.iter().enumerate().for_each(|(idx, gate)| {
                let output = match gate.operation {
                    Operation::Add => running_inputs[gate.left] + running_inputs[gate.right],
                    Operation::Mul => running_inputs[gate.left] * running_inputs[gate.right],
                };

                next_inputs[idx] = output;
            });

            evaluation_layers.push(MultiLinearPolynomial::new(next_inputs.clone()));
            running_inputs = next_inputs;
        });

        let output_polynomial = evaluation_layers.last().unwrap().clone();

        self.layer_evaluations = evaluation_layers;

        output_polynomial
    }

    fn get_bit_idx(
        &self,
        output_idx: usize,
        left_idx: usize,
        right_idx: usize,
        input_bit_repr: usize,
    ) -> usize {
        (((output_idx << input_bit_repr) | left_idx) << input_bit_repr) | right_idx
    }

    fn match_gate_condition(&self, gate: &Gate, condition: &Operation) -> bool {
        match gate.operation {
            Operation::Add => match condition {
                Operation::Add => true,
                _ => false,
            },
            Operation::Mul => match condition {
                Operation::Mul => true,
                _ => false,
            },
        }
    }

    fn get_gate_poly(&self, layer_idx: usize, condition: Operation) -> MultiLinearPolynomial<T> {
        if layer_idx >= self.layers.len() {
            panic!("layer index out of bounds");
        }

        let gates = &self.layers[self.layers.len() - layer_idx - 1];

        let mut output_length = gates.len().next_power_of_two() as usize;

        match output_length {
            1 => output_length = 2 * output_length,
            _ => (),
        }

        let input_lengths_vec: Vec<usize> = gates.iter().fold(vec![], |acc, gate| {
            let mut new_acc = vec![gate.left + 1, gate.right + 1];

            new_acc.extend(&acc);

            new_acc
        });

        let input_bit_length = input_lengths_vec
            .iter()
            .max()
            .unwrap()
            .next_power_of_two()
            .ilog2() as usize;

        let mut evaluation_points: Vec<T> =
            vec![T::from(0); output_length * (1 << (2 * input_bit_length)) as usize];

        gates.iter().enumerate().for_each(|(idx, gate)| {
            if self.match_gate_condition(&gate, &condition) {
                evaluation_points[self.get_bit_idx(idx, gate.left, gate.right, input_bit_length)] =
                    T::from(1);
            }
        });

        MultiLinearPolynomial::new(evaluation_points)
    }

    pub fn get_w_i(&self, layer_idx: usize) -> MultiLinearPolynomial<T> {
        if layer_idx >= self.layer_evaluations.len() {
            panic!("layer index out of bounds");
        }

        self.layer_evaluations[self.layer_evaluations.len() - layer_idx - 1].clone()
    }

    pub fn get_add_i(&self, layer_idx: usize) -> MultiLinearPolynomial<T> {
        self.get_gate_poly(layer_idx, Operation::Add)
    }

    pub fn get_mul_i(&self, layer_idx: usize) -> MultiLinearPolynomial<T> {
        self.get_gate_poly(layer_idx, Operation::Mul)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    fn init_circuit_and_evaluate() -> Circuit<Fq> {
        let mut circuit = Circuit::new(vec![
            vec![
                Gate::new(0, 1, Operation::Add),
                Gate::new(2, 3, Operation::Mul),
            ],
            vec![Gate::new(0, 1, Operation::Add)],
        ]);

        circuit.evaluate_at_input(vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)]);

        circuit
    }

    #[test]
    pub fn arithmetic_gate_test() {
        let circuit = init_circuit_and_evaluate();

        assert_eq!(
            *circuit
                .layer_evaluations
                .last()
                .unwrap()
                .get_evaluation_points(),
            vec![Fq::from(15), Fq::from(0)]
        );
    }

    #[test]
    pub fn test_get_add_i() {
        let circuit = init_circuit_and_evaluate();

        let mut result_vec = vec![Fq::from(0); 32];
        result_vec[1] = Fq::from(1);

        assert_eq!(*circuit.get_add_i(1).get_evaluation_points(), result_vec);
    }

    #[test]
    pub fn test_get_mul_i() {
        let circuit = init_circuit_and_evaluate();

        let mut result_vec = vec![Fq::from(0); 32];
        result_vec[27] = Fq::from(1);

        assert_eq!(*circuit.get_mul_i(1).get_evaluation_points(), result_vec);
    }
}
