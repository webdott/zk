use ark_ff::PrimeField;

struct MultiLinearPolynomial<T: PrimeField> {
    coefficients: Vec<(Vec<T>, T)>,
}

impl<T: PrimeField> MultiLinearPolynomial<T> {
    fn new() -> Self {
        MultiLinearPolynomial {
            coefficients: vec![],
        }
    }

    pub fn partially_evaluate(&self, x: T) -> T {
        todo!()
    }

    pub fn evaluate(&self, x: T) -> T {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
}
