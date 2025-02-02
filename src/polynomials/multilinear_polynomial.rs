use ark_ff::{BigInteger, PrimeField};
use bincode::serialize;
use std::iter;

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
        let serializable_points: Vec<Vec<u8>> = self
            .evaluation_points
            .iter()
            .map(|point| point.into_bigint().to_bytes_le())
            .collect();

        // Serialize the serializable format
        serialize(&serializable_points).expect("Serialization of Multilinear struct failed")
    }

    pub fn evaluation_sum(&self) -> T {
        self.evaluation_points.iter().sum()
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
}
