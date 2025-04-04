use crate::product_polynomial::ProductPolynomial;

use ark_ff::PrimeField;
use std::iter;

#[derive(Debug, Clone)]
pub struct SumPolynomial<T: PrimeField> {
    pub prod_polys: Vec<ProductPolynomial<T>>,
}

impl<T: PrimeField> SumPolynomial<T> {
    pub fn new(prod_polys: Vec<ProductPolynomial<T>>) -> Self {
        let general_poly_length = Self::get_poly_length(&prod_polys);

        prod_polys.iter().for_each(|poly| {
            if poly.length() != general_poly_length {
                panic!("All polynomials must have the same length");
            }
        });

        Self { prod_polys }
    }

    pub fn partial_evaluate(&self, t: &[Option<T>]) -> Self {
        let new_polys = self.prod_polys.iter().map(|poly| poly.partial_evaluate(t));

        Self {
            prod_polys: new_polys.collect(),
        }
    }

    // Evaluate all product polynomials and perform element wise addition
    pub fn evaluate(&self, t: &[Option<T>]) -> T {
        self.prod_polys.iter().map(|poly| poly.evaluate(t)).sum()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.prod_polys
            .iter()
            .flat_map(|poly| poly.to_bytes())
            .collect()
    }

    pub fn reduce(&self) -> Vec<T> {
        // perform element wise product on each multilinear polynomial
        let general_poly_length = self.length();
        let reduced_product_polys: Vec<Vec<T>> =
            self.prod_polys.iter().map(|poly| poly.reduce()).collect();

        let res = iter::repeat(())
            .enumerate()
            .map(|(index, _)| {
                let mut running_idx_sum = T::from(0);

                reduced_product_polys.iter().for_each(|poly| {
                    running_idx_sum += poly[index];
                });

                running_idx_sum
            })
            .take(general_poly_length)
            .collect();

        res
    }

    pub fn get_poly_length(prod_polys: &[ProductPolynomial<T>]) -> usize {
        prod_polys.first().unwrap().length()
    }

    pub fn length(&self) -> usize {
        Self::get_poly_length(&self.prod_polys)
    }

    pub fn number_of_variables(&self) -> u32 {
        self.length().ilog2()
    }

    pub fn degree(&self) -> usize {
        self.prod_polys.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use field_tracker::{print_summary, Ft};

    use crate::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;

    type Fq = Ft!(ark_bn254::Fq);

    fn get_test_prod_polynomial() -> ProductPolynomial<Fq> {
        ProductPolynomial::new(vec![
            MultiLinearPolynomial::new(&vec![Fq::from(2), Fq::from(3), Fq::from(4), Fq::from(5)]),
            MultiLinearPolynomial::new(&vec![Fq::from(2), Fq::from(3), Fq::from(4), Fq::from(5)]),
        ])
    }

    fn get_test_sum_polynomial() -> SumPolynomial<Fq> {
        SumPolynomial::new(vec![get_test_prod_polynomial(), get_test_prod_polynomial()])
    }

    #[test]
    fn test_sum_polynomial_reduce() {
        let test_poly = get_test_sum_polynomial();

        assert_eq!(
            test_poly.reduce(),
            vec![Fq::from(8), Fq::from(18), Fq::from(32), Fq::from(50)]
        );

        print_summary!();
    }

    #[test]
    fn test_sum_polynomial_evaluate() {
        let test_poly = get_test_sum_polynomial();

        assert_eq!(
            test_poly.evaluate(&vec![Some(Fq::from(1)), Some(Fq::from(2))]),
            Fq::from(72)
        );

        print_summary!();
    }

    #[test]
    fn test_sum_polynomial_partial_evaluate() {
        let test_poly = get_test_sum_polynomial();

        assert_eq!(
            test_poly
                .partial_evaluate(&vec![Some(Fq::from(1)), None])
                .prod_polys,
            vec![
                ProductPolynomial::new(vec![
                    MultiLinearPolynomial::new(&vec![Fq::from(4), Fq::from(5)]),
                    MultiLinearPolynomial::new(&vec![Fq::from(4), Fq::from(5)])
                ]),
                ProductPolynomial::new(vec![
                    MultiLinearPolynomial::new(&vec![Fq::from(4), Fq::from(5)]),
                    MultiLinearPolynomial::new(&vec![Fq::from(4), Fq::from(5)])
                ])
            ]
        );

        print_summary!();
    }
}
