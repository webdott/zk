use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use ark_ff::PrimeField;
use std::iter;

#[derive(Debug, Clone)]
pub struct ProductPolynomial<T: PrimeField> {
    polys: Vec<MultiLinearPolynomial<T>>,
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

    pub fn get_multi_polys(&self) -> &Vec<MultiLinearPolynomial<T>> {
        &self.polys
    }

    pub fn partial_evaluate(&self, t: &Vec<Option<T>>) -> Self {
        let new_polys = self.polys.iter().map(|poly| poly.evaluate(t.clone()));

        Self {
            polys: new_polys.collect(),
        }
    }

    pub fn evaluate(&self, t: &[Option<T>]) -> T {
        self.polys
            .iter()
            .map(|poly| {
                poly.evaluate(Vec::from(t))
                    .get_evaluation_points()
                    .first()
                    .unwrap()
                    .clone()
            })
            .product()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.polys.iter().flat_map(|poly| poly.to_bytes()).collect()
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
}

#[cfg(test)]
mod tests {
    #[test]
    fn product_polynomial() {
        todo!()
    }
}
