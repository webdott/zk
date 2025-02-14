use crate::polynomials::product_polynomial::ProductPolynomial;
use ark_ff::PrimeField;
use std::iter;

#[derive(Debug, Clone)]
pub struct SumPolynomial<T: PrimeField> {
    prod_polys: Vec<ProductPolynomial<T>>,
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

    pub fn partial_evaluate(&self, t: &Vec<Option<T>>) -> Self {
        let new_polys = self.prod_polys.iter().map(|poly| poly.partial_evaluate(t));

        Self {
            prod_polys: new_polys.collect(),
        }
    }

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
}

#[cfg(test)]
mod tests {
    #[test]
    fn product_polynomial() {
        todo!()
    }
}
