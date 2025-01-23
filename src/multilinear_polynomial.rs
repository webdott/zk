use ark_ff::PrimeField;
use std::iter;

struct MultiLinearPolynomial<T: PrimeField> {
    evaluation_points: Vec<T>,
}

impl<T: PrimeField> MultiLinearPolynomial<T> {
    fn new(evaluation_points: Vec<T>) -> Self {
        MultiLinearPolynomial { evaluation_points }
    }

    // Given the index where the bit in question is turned off, return flipped index
    fn get_flipped_bit_index(&self, variable_index: usize, init_variable_point: usize) -> usize {
        let power = self.evaluation_points.len().ilog2() - 1 - (variable_index as u32);

        init_variable_point + 2_usize.pow(power)
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

        let mut init_variable_point = 0;
        let power = self.evaluation_points.len() - 1 - variable.0;

        let new_evaluation_points = iter::repeat(())
            .map(|()| {
                let y1 = &self.evaluation_points[init_variable_point];
                let y2 = &self.evaluation_points
                    [self.get_flipped_bit_index(variable.0, init_variable_point)];

                init_variable_point = if (init_variable_point + 1) % 2_usize.pow(power as u32) == 0
                {
                    y2.into_bigint().as_ref()[0] as usize + 1
                } else {
                    init_variable_point + 1
                };

                *y1 + ((*y2 - *y1) * variable.1)
            })
            .take(new_evaluation_points_length)
            .collect();

        MultiLinearPolynomial {
            evaluation_points: new_evaluation_points,
        }
    }

    pub fn evaluate(&self, points: Vec<Option<T>>) -> Self {
        let mut done = 0;

        points
            .iter()
            .enumerate()
            .map(|(idx, point)| (idx, point))
            .fold(
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
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_evaluate() {
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
    pub fn test_evaluate_incomplete_points() {
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
}
