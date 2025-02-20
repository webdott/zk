use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use ark_ff::PrimeField;
use std::iter;

#[derive(Debug, Clone, PartialEq)]
pub struct ProductPolynomial<T: PrimeField> {
    pub polys: Vec<MultiLinearPolynomial<T>>,
}

impl<T: PrimeField> ProductPolynomial<T> {
    pub fn new(polys: Vec<MultiLinearPolynomial<T>>) -> Self {
        let general_poly_length = Self::get_poly_length(&polys);

        polys.iter().for_each(|poly| {
            if poly.get_evaluation_points().len() != general_poly_length {
                panic!("All polynomials must have the same length");
            }
        });

        Self { polys }
    }
    pub fn get_poly_length(polys: &Vec<MultiLinearPolynomial<T>>) -> usize {
        polys.first().unwrap().get_evaluation_points().len()
    }

    pub fn length(&self) -> usize {
        Self::get_poly_length(&self.polys)
    }

    pub fn partial_evaluate(&self, t: &[Option<T>]) -> Self {
        if t.len() != self.length().ilog2() as usize {
            panic!("evaluation points not equal to number of variables");
        }

        let new_polys = self.polys.iter().map(|poly| poly.evaluate(t));

        Self {
            polys: new_polys.collect(),
        }
    }

    pub fn evaluate(&self, t: &[Option<T>]) -> T {
        if t.len() != self.length().ilog2() as usize {
            panic!("evaluation points not equal to number of variables");
        }

        self.polys
            .iter()
            .map(|poly| *poly.evaluate(t).get_evaluation_points().first().unwrap())
            .product()
    }

    pub fn reduce(&self) -> Vec<T> {
        // perform element wise product on each multilinear polynomial
        let general_poly_length = Self::get_poly_length(&self.polys);

        let res = iter::repeat(())
            .enumerate()
            .map(|(index, _)| {
                let mut running_idx_prod = T::from(1);

                self.polys.iter().for_each(|poly| {
                    running_idx_prod *= poly.get_evaluation_points()[index];
                });

                running_idx_prod
            })
            .take(general_poly_length)
            .collect();

        res
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.polys.iter().flat_map(|poly| poly.to_bytes()).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
    use crate::polynomials::product_polynomial::ProductPolynomial;
    use ark_bn254::Fq;

    fn get_test_product_polynomial() -> ProductPolynomial<Fq> {
        ProductPolynomial::new(vec![
            MultiLinearPolynomial::new(&vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)]),
            MultiLinearPolynomial::new(&vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)]),
        ])
    }

    #[test]
    fn test_product_polynomial_reduce() {
        let test_poly = get_test_product_polynomial();

        assert_eq!(
            test_poly.reduce(),
            vec![Fq::from(1), Fq::from(4), Fq::from(9), Fq::from(16)]
        );
    }

    #[test]
    fn test_product_polynomial_evaluate() {
        let test_poly = get_test_product_polynomial();

        assert_eq!(
            test_poly.evaluate(&vec![Some(Fq::from(1)), Some(Fq::from(2))]),
            Fq::from(25)
        );
    }

    #[test]
    fn test_product_polynomial_partial_evaluate() {
        let test_poly = get_test_product_polynomial();

        assert_eq!(
            test_poly
                .partial_evaluate(&vec![Some(Fq::from(1)), None])
                .polys,
            vec![
                MultiLinearPolynomial::new(&vec![Fq::from(3), Fq::from(4)]),
                MultiLinearPolynomial::new(&vec![Fq::from(3), Fq::from(4)])
            ]
        );
    }
}
