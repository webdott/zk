use crate::multilinear::utils::blowup;
use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;

pub struct MultilinearKZGProver<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    _marker_2: PhantomData<P>,
}

impl<T: PrimeField, P: Pairing> MultilinearKZGProver<T, P> {
    fn evaluate_at_tau(
        polynomial: &MultiLinearPolynomial<T>,
        encrypted_lagrange_basis: &[P::G1],
    ) -> P::G1 {
        let polynomial_evals = polynomial.get_evaluation_points();

        if polynomial_evals.len() != encrypted_lagrange_basis.len() {
            panic!("Number of variables of polynomial does not match the number of Taus given!")
        };

        let evaluation_points = (0..encrypted_lagrange_basis.len())
            .map(|i| encrypted_lagrange_basis[i].mul_bigint(polynomial_evals[i].into_bigint()))
            .collect::<Vec<_>>();

        evaluation_points.iter().sum::<P::G1>()
    }

    // Commitment is gotten by doing an element wise multiplication between encrypted lagrange basis and the multilinear polynomial
    pub fn generate_commitment(
        f: &MultiLinearPolynomial<T>,
        encrypted_lagrange_basis: &[P::G1],
    ) -> P::G1 {
        Self::evaluate_at_tau(f, encrypted_lagrange_basis)
    }

    pub fn generate_proof(
        openings: &[T],
        encrypted_lagrange_basis: &[P::G1],
        polynomial: &MultiLinearPolynomial<T>,
    ) -> (T, Vec<P::G1>) {
        // divide the polynomial by each opening as a factor
        // e.g. if the roots are a = 6, b = 7, c = 0; we divide the polynomial by a-6, remainder by b - 7 and lastly, c - 0;
        // But in actual fact, we are evaluating the polynomial at the variable points.

        let opening_points = openings.iter().map(|val| Some(*val)).collect::<Vec<_>>();
        let v_poly = polynomial.evaluate(&opening_points);

        let mut quotients = Vec::with_capacity(openings.len());

        let f_minus_v = polynomial.minus(v_poly.get_evaluation_points().first().unwrap());
        let mut dividend = f_minus_v;

        for (idx, opening) in openings.iter().enumerate() {
            let (mut quotient, remainder) = dividend.compute_quotient_remainder(opening, idx);

            dividend = remainder;
            quotient = blowup(openings.len(), idx, &quotient);
            let quotient_at_tau: P::G1 = Self::evaluate_at_tau(
                &MultiLinearPolynomial::new(&quotient),
                encrypted_lagrange_basis,
            );

            quotients.push(quotient_at_tau);
        }

        (*v_poly.get_evaluation_points().first().unwrap(), quotients)
    }
}

#[cfg(test)]
mod tests {
    use crate::multilinear::utils::{
        encrypt_lagrange_basis, generate_lagrange_basis_for_n_variables,
    };

    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Affine};
    use ark_ec::AffineRepr;

    #[test]
    pub fn test_generate_commitment() {
        let encrypted_lagrange_basis = encrypt_lagrange_basis::<Fr, Bls12_381>(
            &generate_lagrange_basis_for_n_variables(3, &[Fr::from(5), Fr::from(2), Fr::from(3)]),
        );

        let commitment = MultilinearKZGProver::<Fr, Bls12_381>::generate_commitment(
            &MultiLinearPolynomial::new(&vec![
                Fr::from(0),
                Fr::from(4),
                Fr::from(0),
                Fr::from(4),
                Fr::from(0),
                Fr::from(4),
                Fr::from(3),
                Fr::from(7),
            ]),
            &encrypted_lagrange_basis,
        );

        assert_eq!(
            commitment,
            G1Affine::generator().mul_bigint(Fr::from(42).into_bigint())
        )
    }
}
