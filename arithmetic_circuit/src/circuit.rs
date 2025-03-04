use crate::gate::{Gate, Operation};

use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;

use ark_ff::PrimeField;
use std::cmp::max;
use std::marker::PhantomData;

pub struct Circuit<T: PrimeField> {
    _marker: PhantomData<T>,
    layers: Vec<Vec<Gate>>,
}

impl<T: PrimeField> Circuit<T> {
    pub fn new(layers: Vec<Vec<Gate>>) -> Self {
        Self {
            layers,
            _marker: PhantomData,
        }
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

    // This takes in a set of inputs and for each layer of gates we have, calculate the next set of inputs
    // The set of inputs are stored as evaluation layers for easy retrieval
    pub fn evaluate_at_input(&mut self, inputs: Vec<T>) -> Vec<MultiLinearPolynomial<T>> {
        let mut evaluation_layers = vec![MultiLinearPolynomial::new(&inputs)];
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

            evaluation_layers.push(MultiLinearPolynomial::new(&next_inputs));
            running_inputs = next_inputs;
        });

        evaluation_layers
    }

    // This helps us to get the index at which a gate is present (turned on)
    // Say we have a gate that has an output at index 00, left input at 10 and right input at 11
    // The index for that in the gate poly would be 001011
    // To achieve this, we use bit manipulation - combining the left shift and OR operations.
    // Left shift to accommodate for next index to add and OR operator to add the index.
    fn get_bit_idx(
        &self,
        output_idx: usize,
        left_idx: usize,
        right_idx: usize,
        input_bit_repr: usize,
    ) -> usize {
        (((output_idx << input_bit_repr) | left_idx) << input_bit_repr) | right_idx
    }

    // This gets the gate polynomial at an index represented in multilinear form
    // For each gate have an output index, two input indexes for the two inputs
    // In this case, the output is basically the index of the gate since they are in a vec
    // The evaluation points would basically be 2^(all bits used to represent output, and the two indexes).
    // i.e if we have the gate at output index 10, left input index at 00 and right index at 01:
    // In total, there are 6 bits (100001) in total used to represent this gate poly which is 2^6 evaluation points.
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
                // set the index where gate is present to 1.
                evaluation_points[self.get_bit_idx(idx, gate.left, gate.right, input_bit_length)] =
                    T::from(1);
            }
        });

        MultiLinearPolynomial::new(&evaluation_points)
    }

    // After evaluation of the circuit, we can just get the polynomial of each w_layer
    pub fn get_w_i(
        &self,
        layer_idx: usize,
        layer_evaluations: &[MultiLinearPolynomial<T>],
    ) -> MultiLinearPolynomial<T> {
        if layer_idx >= layer_evaluations.len() {
            panic!("layer index out of bounds");
        }

        MultiLinearPolynomial::new(
            layer_evaluations[layer_evaluations.len() - layer_idx - 1].get_evaluation_points(),
        )
    }

    pub fn get_add_i(&self, layer_idx: usize) -> MultiLinearPolynomial<T> {
        self.get_gate_poly(layer_idx, Operation::Add)
    }

    pub fn get_mul_i(&self, layer_idx: usize) -> MultiLinearPolynomial<T> {
        self.get_gate_poly(layer_idx, Operation::Mul)
    }

    // Calculate how many layers we have in the circuit
    pub fn get_layer_count(&self) -> usize {
        self.layers.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    fn init_circuit_and_evaluate() -> (Vec<MultiLinearPolynomial<Fq>>, Circuit<Fq>) {
        let mut circuit = Circuit::new(vec![
            vec![
                Gate::new(0, 1, Operation::Add),
                Gate::new(2, 3, Operation::Mul),
            ],
            vec![Gate::new(0, 1, Operation::Add)],
        ]);

        (
            circuit.evaluate_at_input(vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)]),
            circuit,
        )
    }

    #[test]
    pub fn arithmetic_gate_test() {
        let (circuit_evaluations, _) = init_circuit_and_evaluate();

        assert_eq!(
            *circuit_evaluations.last().unwrap().get_evaluation_points(),
            vec![Fq::from(15), Fq::from(0)]
        );
    }

    #[test]
    pub fn test_get_add_i() {
        let (_, circuit) = init_circuit_and_evaluate();

        let mut result_vec = vec![Fq::from(0); 32];
        result_vec[1] = Fq::from(1);

        assert_eq!(*circuit.get_add_i(1).get_evaluation_points(), result_vec);
    }

    #[test]
    pub fn test_get_mul_i() {
        let (_, circuit) = init_circuit_and_evaluate();

        let mut result_vec = vec![Fq::from(0); 32];
        result_vec[27] = Fq::from(1);

        assert_eq!(*circuit.get_mul_i(1).get_evaluation_points(), result_vec);
    }
}
