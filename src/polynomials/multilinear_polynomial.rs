use ark_ff::{BigInteger, PrimeField};
use std::iter;

use crate::concepts_protocols::arithmetic_gate::gate::Operation;

#[derive(Debug, Clone)]
pub struct MultiLinearPolynomial<T: PrimeField> {
    evaluation_points: Vec<T>,
}

impl<T: PrimeField> MultiLinearPolynomial<T> {
    pub fn new(evaluation_points: Vec<T>) -> Self {
        if !evaluation_points.len().is_power_of_two() {
            panic!(
                "Invalid Multilinear Polynomial: evaluation points length is not a power of two"
            );
        }

        Self { evaluation_points }
    }

    pub fn get_evaluation_points(&self) -> &Vec<T> {
        &self.evaluation_points
    }

    pub fn number_of_variables(&self) -> u32 {
        self.evaluation_points.len().ilog2()
    }

    // Given the index where the bit in question is turned off, return flipped index
    fn get_flipped_bit_with_bitwise_or(
        &self,
        index_to_flip: usize,
        number_to_flip: usize,
    ) -> usize {
        let power = self.number_of_variables() - 1 - (index_to_flip as u32);

        number_to_flip | (1 << power)
    }

    pub fn partially_evaluate(&self, variable: (usize, T)) -> Self {
        // interpolate + evaluate
        // y1 + r(y2 - y1)
        // where:
        // y1 -> first point
        // y2 -> second point
        // r -> evaluation point

        let evaluation_length = self.evaluation_points.len();
        let new_evaluation_points_length = evaluation_length / 2;

        let mut y1_index = 0;
        let target = 1 << (self.number_of_variables() as usize) - 1 - variable.0;

        let new_evaluation_points = iter::repeat(())
            .map(|()| {
                let y1 = &self.evaluation_points[y1_index];
                let y2_index = self.get_flipped_bit_with_bitwise_or(variable.0, y1_index);
                let y2 = &self.evaluation_points[y2_index];

                y1_index = if (y1_index + 1) % target == 0 {
                    y2_index + 1
                } else {
                    y1_index + 1
                };

                *y1 + ((*y2 - *y1) * variable.1)
            })
            .take(new_evaluation_points_length)
            .collect();

        Self::new(new_evaluation_points)
    }

    pub fn evaluate(&self, points: Vec<Option<T>>) -> Self {
        if points.len() != self.number_of_variables() as usize {
            panic!("points length does not match number of variables");
        }

        let mut done = 0;

        points.iter().enumerate().fold(
            MultiLinearPolynomial::new(self.evaluation_points.clone()),
            |acc, (idx, point)| {
                let mlp = match point {
                    Some(_) => {
                        let new_acc = acc.partially_evaluate((idx - done, point.unwrap()));
                        done += 1;

                        new_acc
                    }
                    None => acc,
                };

                mlp
            },
        )
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        // Convert evaluation points to a serializable format (e.g., bytes)
        let serializable_points: Vec<u8> = self
            .evaluation_points
            .iter()
            .flat_map(|point| point.into_bigint().to_bytes_le())
            .collect();

        serializable_points
    }

    pub fn evaluation_sum(&self) -> T {
        self.evaluation_points.iter().sum()
    }

    fn operate_w(
        w_b: &MultiLinearPolynomial<T>,
        w_c: &MultiLinearPolynomial<T>,
        operation: Operation,
    ) -> MultiLinearPolynomial<T> {
        let (evals_b, evals_c) = (w_b.get_evaluation_points(), w_c.get_evaluation_points());
        let (w_b_len, w_c_len) = (evals_b.len(), evals_c.len());
        let result_evaluation_length = w_b_len * w_c_len;
        let mut result_eval_points = vec![T::from(0); result_evaluation_length];

        vec![0; result_evaluation_length]
            .iter()
            .enumerate()
            .for_each(|(i, _)| {
                let (idx_b, idx_c) = ((i / w_b_len), i % w_c_len);

                match operation {
                    Operation::Add => {
                        result_eval_points[i] = evals_b[idx_b] + evals_c[idx_c];
                    }
                    Operation::Mul => result_eval_points[i] = evals_b[idx_b] * evals_c[idx_c],
                }
            });

        Self::new(result_eval_points)
    }

    pub fn w_add(
        w_b: &MultiLinearPolynomial<T>,
        w_c: &MultiLinearPolynomial<T>,
    ) -> MultiLinearPolynomial<T> {
        Self::operate_w(w_b, w_c, Operation::Add)
    }

    pub fn w_mul(&self, other: &MultiLinearPolynomial<T>) -> MultiLinearPolynomial<T> {
        Self::operate_w(other, self, Operation::Mul)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_evaluate_4_variables() {
        // 3ac + 4bd + 5ab -> where a = 4, b = 2, c = 6, d = 1
        let mlp = MultiLinearPolynomial::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(9),
            Fq::from(8),
            Fq::from(12),
        ]);

        assert_eq!(
            mlp.evaluate(vec![
                Some(Fq::from(4)),
                Some(Fq::from(2)),
                Some(Fq::from(6)),
                Some(Fq::from(1))
            ])
            .evaluation_points,
            vec![Fq::from(120)],
        );
    }

    #[test]
    pub fn test_partially_evaluate_4_variables_incomplete_points() {
        // 3ac + 4bd + 5ab -> where a = 4
        let mlp = MultiLinearPolynomial::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(9),
            Fq::from(8),
            Fq::from(12),
        ]);

        assert_eq!(
            mlp.evaluate(vec![Some(Fq::from(4)), None, None, None])
                .evaluation_points,
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(12),
                Fq::from(12),
                Fq::from(20),
                Fq::from(24),
                Fq::from(32),
                Fq::from(36)
            ],
        );
    }

    #[test]
    pub fn test_partially_evaluate_3_variables_incomplete_points() {
        //2ab + 3bc -> where c = 3
        let mlp = MultiLinearPolynomial::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ]);

        assert_eq!(
            mlp.evaluate(vec![None, None, Some(Fq::from(3))])
                .evaluation_points,
            vec![Fq::from(0), Fq::from(9), Fq::from(0), Fq::from(11)],
        );
    }

    #[test]
    pub fn test_operation_ws() {
        let (w_b, w_c) = (
            MultiLinearPolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]),
            MultiLinearPolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]),
        );

        assert_eq!(
            *MultiLinearPolynomial::w_add(&w_b, &w_c).get_evaluation_points(),
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(2),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(2),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(2),
                Fq::from(3),
                Fq::from(3),
                Fq::from(3),
                Fq::from(5)
            ]
        );

        assert_eq!(
            *MultiLinearPolynomial::w_mul(&w_b, &w_c).get_evaluation_points(),
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(6)
            ]
        );
    }
}
