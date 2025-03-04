use ark_ff::{BigInteger, PrimeField};
use std::ops::Add;

#[derive(Debug)]
enum Operation {
    Add,
    Mul,
}

pub enum BlowUpDirection {
    Left,
    Right,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MultiLinearPolynomial<T: PrimeField> {
    evaluation_points: Vec<T>,
}

impl<T: PrimeField> MultiLinearPolynomial<T> {
    pub fn new(evaluation_points: &Vec<T>) -> Self {
        if !evaluation_points.len().is_power_of_two() {
            panic!(
                "Invalid Multilinear Polynomial: evaluation points length is not a power of two"
            );
        }

        Self {
            evaluation_points: evaluation_points.to_vec(),
        }
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

    fn get_y1_y2_indexes(&self, variable_idx: usize) -> Vec<(usize, usize)> {
        // This gets the different pairing indexes for y1 and y2.
        // Starting from 0, we can always get the next bit at which that bit is turned on
        // We basically get all the indexes where the position of the variable is turned off, set that to y1,
        // We also find where it's turned on and set that to y2
        // E.g. you have variables a, b, c; Boolean hypercube values are below. For variable b:
        // 0 0 0, -> y1_0
        // 0 0 1, -> y1_1
        // 0 1 0, -> y2_0
        // 0 1 1, -> y2_1
        // 1 0 0, -> y1_2
        // 1 0 1, -> y1_3
        // 1 1 0, -> y2_2
        // 1 1 1, -> y2_3
        // From pattern, you can see that we skip every 2^power_index of variable.

        let mut y1_index = 0;
        let half_length = self.evaluation_points.len() / 2;
        let mut y1_y2_indexes = Vec::with_capacity(half_length);
        let target = 1 << ((self.number_of_variables() as usize) - 1 - variable_idx);

        (0..half_length).for_each(|_| {
            let y2_index = self.get_flipped_bit_with_bitwise_or(variable_idx, y1_index);

            y1_y2_indexes.push((y1_index, y2_index));

            // Check if we should skip for next y1 or just continue to the next one
            y1_index = if (y1_index + 1) % target == 0 {
                y2_index + 1
            } else {
                y1_index + 1
            };
        });

        y1_y2_indexes
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

        // This performs tensor addition or multiplication between two polynomials of different variables
        // This variant uses tensor addition and multiplication:
        // W(b) * W(c) =>
        // [2, 2] * [3, 2] => [2 * 3, 2 * 2, 2 * 3, 2 * 2] => [6, 4, 6, 4]
        (0..result_evaluation_length)
            .enumerate()
            .for_each(|(i, _)| {
                let (idx_b, idx_c) = (i / w_b_len, i % w_c_len);

                match operation {
                    Operation::Add => result_eval_points[i] = evals_b[idx_b] + evals_c[idx_c],
                    Operation::Mul => result_eval_points[i] = evals_b[idx_b] * evals_c[idx_c],
                }
            });

        Self::new(&result_eval_points)
    }

    // This repeats the evaluation points as a whole twice.
    // Equivalent to adding 0 as the coefficient of missing variables at the back
    // E.g. 4b => 0a + 4b
    fn blowup_left(evaluation_points: &[T]) -> Vec<T> {
        let mut blown_up_points = Vec::from(evaluation_points);
        blown_up_points.extend_from_slice(evaluation_points);

        blown_up_points
    }

    // This repeats each evaluation point twice side-by-side
    // Equivalent to adding 0 as the coefficient of missing variables in front
    // E.g. 4a => 4a + 0b
    fn blowup_right(evaluation_points: &[T]) -> Vec<T> {
        let mut blown_up_points = Vec::with_capacity(evaluation_points.len() * 2);

        evaluation_points.iter().for_each(|e| {
            blown_up_points.push(*e);
            blown_up_points.push(*e);
        });

        blown_up_points
    }

    pub fn blow_up_n_times(dir: BlowUpDirection, n: usize, evaluation_points: &[T]) -> Vec<T> {
        let mut running_evaluation_points = Vec::from(evaluation_points);

        for _ in 0..n {
            running_evaluation_points = match dir {
                BlowUpDirection::Left => Self::blowup_left(&running_evaluation_points),
                BlowUpDirection::Right => Self::blowup_right(&running_evaluation_points),
            }
        }

        running_evaluation_points
    }

    pub fn scalar_mul(&self, scalar: T) -> Self {
        Self::new(&self.evaluation_points.iter().map(|e| *e * scalar).collect())
    }

    pub fn get_evaluation_points(&self) -> &Vec<T> {
        &self.evaluation_points
    }

    pub fn number_of_variables(&self) -> u32 {
        self.evaluation_points.len().ilog2()
    }

    pub fn partially_evaluate(&self, variable: (usize, T)) -> Self {
        // interpolate + evaluate
        // y1 + r(y2 - y1)
        // where:
        // y1 -> first point
        // y2 -> second point
        // r -> evaluation point
        let new_evaluation_points_length = self.evaluation_points.len() / 2;
        let y1_y2_indexes = self.get_y1_y2_indexes(variable.0);

        // Given the various pairing indexes for y1 and y2, carry out formula
        let new_evaluation_points = y1_y2_indexes
            .iter()
            .map(|(y1_index, y2_index)| {
                let (y1, y2) = (
                    self.evaluation_points[*y1_index],
                    self.evaluation_points[*y2_index],
                );

                y1 + ((y2 - y1) * variable.1)
            })
            .take(new_evaluation_points_length)
            .collect();

        Self::new(&new_evaluation_points)
    }

    pub fn evaluate(&self, points: &[Option<T>]) -> Self {
        if points.len() != self.number_of_variables() as usize {
            panic!("points length does not match number of variables");
        }

        let mut done = 0;

        points.iter().enumerate().fold(
            MultiLinearPolynomial::new(&self.evaluation_points),
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

    // Adds two polynomials of same variables together
    pub fn _add(&self, other: &MultiLinearPolynomial<T>) -> Self {
        if self.number_of_variables() != other.number_of_variables() {
            panic!("Polynomial must have the same length");
        };

        let mut new_evals = vec![T::from(0); other.evaluation_points.len()];

        (0..self.evaluation_points.len()).for_each(|idx| {
            new_evals[idx] +=
                self.get_evaluation_points()[idx] + other.get_evaluation_points()[idx];
        });

        Self::new(&new_evals)
    }

    // Performs F(x) - V operation
    pub fn minus(&self, other: &T) -> Self {
        let new_evaluation_points = self
            .evaluation_points
            .iter()
            .map(|val| *val - *other)
            .collect::<Vec<_>>();

        Self::new(&new_evaluation_points)
    }

    // Dividing a polynomial at a variable point gives you the *quotient* and *remainder*
    // The quotient can be gotten by finding the two points at which the variable switches (from 0 to 1)
    //    - Then subtract the evaluation points at the indexes (y2 - y1)

    // The remainder can be gotten by partially evaluating the polynomial at the variables r.
    pub fn compute_quotient_remainder(&self, divisor: &T, variable_index: usize) -> (Vec<T>, Self) {
        let y1_y2_indexes = self.get_y1_y2_indexes(variable_index);

        let remainder = self.partially_evaluate((variable_index, *divisor));
        let quotient = y1_y2_indexes
            .iter()
            .map(|(y1_index, y2_index)| {
                let (y1, y2) = (
                    self.evaluation_points[*y1_index],
                    self.evaluation_points[*y2_index],
                );

                y2 - y1
            })
            .collect::<Vec<_>>();

        (quotient, remainder)
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

impl<T: PrimeField> Add for MultiLinearPolynomial<T> {
    type Output = Self;

    fn add(self, other: MultiLinearPolynomial<T>) -> Self {
        self._add(&other)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    fn get_test_polynomial() -> MultiLinearPolynomial<Fq> {
        MultiLinearPolynomial::new(&vec![
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
        ])
    }

    fn get_test_polynomial_2() -> MultiLinearPolynomial<Fq> {
        MultiLinearPolynomial::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ])
    }

    #[test]
    pub fn test_evaluate_4_variables() {
        // 3ac + 4bd + 5ab -> where a = 4, b = 2, c = 6, d = 1
        let mlp = get_test_polynomial();

        assert_eq!(
            mlp.evaluate(&vec![
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
        let mlp = get_test_polynomial();

        assert_eq!(
            mlp.evaluate(&vec![Some(Fq::from(4)), None, None, None])
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
        let mlp = get_test_polynomial_2();

        assert_eq!(
            mlp.evaluate(&vec![None, None, Some(Fq::from(3))])
                .evaluation_points,
            vec![Fq::from(0), Fq::from(9), Fq::from(0), Fq::from(11)],
        );
    }

    #[test]
    pub fn test_operation_ws() {
        let (w_b, w_c) = (
            MultiLinearPolynomial::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]),
            MultiLinearPolynomial::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]),
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

    #[test]
    pub fn test_blowup() {
        let init_evaluation_points = &[Fq::from(0), Fq::from(3)];
        let blown_up_left_values = MultiLinearPolynomial::blow_up_n_times(
            BlowUpDirection::Left,
            1,
            init_evaluation_points,
        );

        assert_eq!(
            blown_up_left_values,
            vec![Fq::from(0), Fq::from(3), Fq::from(0), Fq::from(3)]
        );

        assert_eq!(
            MultiLinearPolynomial::blow_up_n_times(
                BlowUpDirection::Right,
                1,
                &blown_up_left_values
            ),
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(3),
                Fq::from(3),
                Fq::from(0),
                Fq::from(0),
                Fq::from(3),
                Fq::from(3),
            ]
        )
    }
}
