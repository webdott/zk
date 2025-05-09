use ark_ff::{BigInteger, PrimeField};
use field_tracker::{end_tscope, start_tscope};
use std::ops::{Add, Mul};
use std::{cmp, mem};

#[derive(Debug)]
pub struct UnivariatePolynomial<T: PrimeField> {
    pub coefficients: Vec<T>,
}

impl<T: PrimeField> UnivariatePolynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        UnivariatePolynomial { coefficients }
    }

    // Given a point, evaluate the result of the polynomial at that point
    // x^2 + 5x + 2 (@ x = 2) => (2 * x^0) + (5 * x ) + (1 * x * x)
    // From this, we can see that rather than raising x to the power each time,
    // we could keep a running product to multiply with the polynomials coefficients
    pub fn evaluate(&self, x: T) -> T {
        start_tscope!("Univariate Polynomial Evaluate");

        let mut result: T = T::from(0);
        let mut running_x: T = T::from(1);

        for i in 0..self.coefficients.len() {
            result += self.coefficients[i] * (running_x);
            running_x *= x;
        }

        end_tscope!();

        result
    }

    // Get evaluation of the polynomial over the boolean hypercube and return sum
    pub fn evaluate_sum_over_boolean_hypercube(&self) -> T {
        start_tscope!("Univariate Polynomial Sum Over Boolean HC");

        let sum = self.evaluate(T::from(0)) + self.evaluate(T::from(1));

        end_tscope!();

        sum
    }

    // Given a specific list of points, find the original polynomial
    // To do this, we use perform Lagrange interpolation on the points:
    //           (x - x1)(x - x2)...(x - xn)                  (x - x0)(x - x2)...(x - xn)                   (x - x0)(x - x1)...(x - xn-1)
    // f(x) =   ------------------------------   * y0   +  ------------------------------   * y1 ..... +  ------------------------------  * yn
    //          (x0 - x1)(x0 - x2)...(x0 - xn)              (x1 - x0)(x1 - x2)...(x1 - xn)                (xn - x0)(xn - x1)...(xn - xn-1)
    pub fn interpolate(x_points: &[T], y_points: &[T]) -> Self {
        start_tscope!("Univariate Interpolate");

        let n = x_points.len();

        let mut res = UnivariatePolynomial {
            coefficients: vec![T::from(0); n],
        };

        for i in 0..n {
            let mut denominator: T = T::from(1);

            // numerator is a multiplication of polynomials together e.g (x - x1)(x - x2)...(x - xn)
            let mut numerator = UnivariatePolynomial {
                coefficients: vec![T::from(1)],
            };

            for j in 0..n {
                if i == j {
                    continue;
                }

                let int_poly = UnivariatePolynomial {
                    coefficients: vec![-T::from(x_points[j]), T::from(1)],
                };

                denominator *= T::from(x_points[i]) - T::from(x_points[j]);
                numerator = numerator * int_poly
            }

            res = res + (numerator.scalar_mul(y_points[i] / denominator));
        }

        end_tscope!();

        res
    }

    // perform scalar mul between number and polynomial. Alternatively, you could represent a constant number as a polynomial i.e
    // UnivariatePolynomial {
    //      coefficients: [1]
    // } which represents the number 1.
    pub fn scalar_mul(&self, num: T) -> Self {
        start_tscope!("Univariate Scalar Mul");

        let n = self.coefficients.len();
        let mut res = vec![T::from(0); n];

        for i in 0..n {
            res[i] = self.coefficients[i] * num;
        }

        end_tscope!();

        UnivariatePolynomial { coefficients: res }
    }

    // Multiply polynomials together
    // You get a polynomial with a degree of the highest degrees in each polynomial multiplied together
    pub fn _mul(&self, p2: &Self) -> Self {
        start_tscope!("Univariate Mul");

        let len_1 = self.coefficients.len();
        let len_2 = p2.coefficients.len();

        let max_len = cmp::max(len_1, len_2);

        let mut greater_coef = &self.coefficients;
        let mut lesser_coef = &p2.coefficients;

        if len_2 > len_1 {
            mem::swap(&mut greater_coef, &mut lesser_coef);
        }

        let mut coefs: Vec<T> = vec![T::from(0); len_2 + len_1 - 1];

        // (2 + x) * (3 + x^2) =>
        // [2, 1] *  [3, 0, 1] =>
        // [6, 0 + 3, 2 + 0, 1] =>
        // [6, 3, 3, 1] => 6 + 3x + 3x^2 + x^3
        for i in 0..max_len {
            let mut idx = i;

            for j in 0..lesser_coef.len() {
                coefs[idx] += greater_coef[i] * lesser_coef[j];

                idx += 1;
            }
        }

        end_tscope!();

        UnivariatePolynomial {
            coefficients: coefs,
        }
    }

    // Add polynomials together
    // You get a polynomial with a degree of the highest degree in any of the polynomial
    pub fn _add(&self, p2: &Self) -> Self {
        start_tscope!("Univariate Add");

        let len_1 = self.coefficients.len();
        let len_2 = p2.coefficients.len();

        let max_len = cmp::min(len_1, len_2);

        let mut coefs = vec![T::from(0); max_len];

        // (2 + x) + (3 + x^2) =>
        // [2, 1]  +  [3, 0, 1] =>
        // [3 + 2, 1 + 0, 1] =>
        // [5, 1, 1] => 5 + x + x^2
        for i in 0..max_len {
            if i < len_1 {
                coefs[i] += self.coefficients[i];
            }

            if i < len_2 {
                coefs[i] += p2.coefficients[i];
            }
        }

        end_tscope!();

        UnivariatePolynomial {
            coefficients: coefs,
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        // Convert evaluation points to a serializable format (e.g., bytes)
        let serializable_points: Vec<u8> = self
            .coefficients
            .iter()
            .flat_map(|point| point.into_bigint().to_bytes_le())
            .collect();

        serializable_points
    }
}

impl<T: PrimeField> Add for UnivariatePolynomial<T> {
    type Output = Self;

    fn add(self, other: UnivariatePolynomial<T>) -> Self {
        self._add(&other)
    }
}

impl<T: PrimeField> Mul for UnivariatePolynomial<T> {
    type Output = Self;

    fn mul(self, p2: UnivariatePolynomial<T>) -> Self {
        self._mul(&p2)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use field_tracker::{print_summary, Ft};

    type Fq = Ft!(ark_bn254::Fq);

    #[test]
    pub fn test_evaluate() {
        let poly = UnivariatePolynomial::new(vec![Fq::from(20), Fq::from(10), Fq::from(3)]);

        assert_eq!(poly.evaluate(Fq::from(2)), Fq::from(52));

        print_summary!();
    }

    #[test]
    pub fn test_mul() {
        let poly1 = UnivariatePolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 =
            UnivariatePolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(4)]);

        assert_eq!(
            (poly1 * poly2).coefficients,
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(0),
                Fq::from(8)
            ]
        );
    }

    #[test]
    pub fn test_fibonacci_range() {
        let poly = UnivariatePolynomial::interpolate(
            &vec![
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
            &vec![
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

        assert_eq!(
            poly.evaluate(Fq::from(5)),
            poly.evaluate(Fq::from(4)) + poly.evaluate(Fq::from(3))
        );

        print_summary!();
    }

    #[test]
    pub fn test_interpolate() {
        let poly1 = UnivariatePolynomial::interpolate(
            &vec![Fq::from(0), Fq::from(1), Fq::from(2)],
            &vec![Fq::from(8), Fq::from(10), Fq::from(16)],
        );

        // f(x) = 2x
        // [(0, 0), (1,2)]
        let poly2 = UnivariatePolynomial::interpolate(
            &vec![Fq::from(0), Fq::from(1)],
            &vec![Fq::from(0), Fq::from(2)],
        );

        // f(x) = 2x
        // [(2, 4), (4,8)]
        let poly3 = UnivariatePolynomial::interpolate(
            &vec![Fq::from(2), Fq::from(4)],
            &vec![Fq::from(4), Fq::from(8)],
        );

        assert_eq!(
            poly1.coefficients,
            vec![Fq::from(8), Fq::from(0), Fq::from(2)]
        );
        assert_eq!(poly2.coefficients, vec![Fq::from(0), Fq::from(2)]);
        assert_eq!(poly3.coefficients, vec![Fq::from(0), Fq::from(2)]);

        print_summary!();
    }
}
