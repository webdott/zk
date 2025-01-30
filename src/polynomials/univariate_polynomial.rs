use ark_ff::PrimeField;
use std::{cmp, mem};

#[derive(Debug)]
pub struct UnivariatePolynomial<T: PrimeField> {
    pub coefficients: Vec<T>,
}

impl<T: PrimeField> UnivariatePolynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        UnivariatePolynomial { coefficients }
    }

    // given a point, evaluate the result of the polynomial at that point
    pub fn evaluate(&self, x: T) -> T {
        let mut result: T = T::from(0);

        for i in 0..self.coefficients.len() {
            result += self.coefficients[i] * (x.pow([i as u64]));
        }

        result
    }

    // Given a specific list of points, find the original polynomial
    pub fn interpolate(x_points: Vec<T>, y_points: Vec<T>) -> Self {
        let n = x_points.len();

        let mut res = UnivariatePolynomial {
            coefficients: vec![T::from(0); n],
        };

        for i in 0..n {
            let mut denominator: T = T::from(1);
            let mut running_poly = UnivariatePolynomial {
                coefficients: vec![T::from(1)],
            };

            for j in 0..n {
                if i == j {
                    continue;
                }

                let int_poly = UnivariatePolynomial {
                    coefficients: vec![-(T::from(x_points[j])) / T::from(1), T::from(1)],
                };

                denominator *= T::from(x_points[i]) - T::from(x_points[j]);
                running_poly = running_poly.mul(&int_poly)
            }

            res = res.add(&running_poly.scalar_mul(y_points[i] / denominator));
        }

        res
    }

    // perform scalar mul between number and polynomial. Alternatively, you could represent a constant number as a polynomial i.e
    // UnivariatePolynomial {
    //      coefficients: [1]
    // } which represents the number 1.
    pub fn scalar_mul(&self, num: T) -> Self {
        let n = self.coefficients.len();
        let mut res = vec![T::from(0); n];

        for i in 0..n {
            res[i] = self.coefficients[i] * num;
        }

        UnivariatePolynomial { coefficients: res }
    }

    // Multiply polynomials together
    pub fn mul(&self, p2: &Self) -> Self {
        let len_1 = self.coefficients.len();
        let len_2 = p2.coefficients.len();

        let max_len = cmp::max(len_1, len_2);

        let mut greater_coef = &self.coefficients;
        let mut lesser_coef = &p2.coefficients;

        if len_2 > len_1 {
            mem::swap(&mut greater_coef, &mut lesser_coef);
        }

        let mut coefs: Vec<T> = vec![T::from(0); max_len + 1];

        for i in 0..max_len {
            let mut idx = i;

            for j in 0..lesser_coef.len() {
                coefs[idx] += greater_coef[i] * lesser_coef[j];

                idx += 1;
            }
        }

        UnivariatePolynomial {
            coefficients: coefs,
        }
    }

    // add polynomials together
    pub fn add(&self, p2: &Self) -> Self {
        let len_1 = self.coefficients.len();
        let len_2 = p2.coefficients.len();

        let max_len = cmp::min(len_1, len_2);

        let mut coefs = vec![T::from(0); max_len];

        for i in 0..max_len {
            if i < len_1 {
                coefs[i] += self.coefficients[i];
            }

            if i < len_2 {
                coefs[i] += p2.coefficients[i];
            }
        }

        UnivariatePolynomial {
            coefficients: coefs,
        }
    }
}

#[cfg(test)]
mod test {
    use super::UnivariatePolynomial;
    use ark_bn254::Fq;
    use rand::Rng;

    #[test]
    pub fn test_interpolate() {
        let poly = UnivariatePolynomial::interpolate(
            vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
                Fq::from(7),
                Fq::from(8),
            ],
            vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(5),
                Fq::from(8),
                Fq::from(13),
                Fq::from(21),
            ],
        );

        let mut random_x = rand::thread_rng();
        let random_u64: u64 = random_x.gen_range(2..100);

        assert_eq!(
            poly.evaluate(Fq::from(random_u64)),
            poly.evaluate(Fq::from(random_u64 - 1)) + poly.evaluate(Fq::from(random_u64 - 2))
        );
    }

    #[test]
    pub fn test_interpolate_2() {
        let poly1 = UnivariatePolynomial::interpolate(
            vec![Fq::from(0), Fq::from(1), Fq::from(2)],
            vec![Fq::from(8), Fq::from(10), Fq::from(16)],
        );

        // f(x) = 2x
        // [(0, 0), (1,2)]
        let poly2 = UnivariatePolynomial::interpolate(
            vec![Fq::from(0), Fq::from(1)],
            vec![Fq::from(0), Fq::from(2)],
        );

        // f(x) = 2x
        // [(2, 4), (4,8)]
        let poly3 = UnivariatePolynomial::interpolate(
            vec![Fq::from(2), Fq::from(4)],
            vec![Fq::from(4), Fq::from(8)],
        );

        assert_eq!(
            poly1.coefficients,
            vec![Fq::from(8), Fq::from(0), Fq::from(2)]
        );
        assert_eq!(poly2.coefficients, vec![Fq::from(0), Fq::from(2)]);
        assert_eq!(poly3.coefficients, vec![Fq::from(0), Fq::from(2)]);
    }
}
