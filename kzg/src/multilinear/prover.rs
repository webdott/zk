use polynomials::multilinear_polynomial::evaluation_form::{
    BlowUpDirection, MultiLinearPolynomial,
};
use std::cmp::max;

use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;
use std::marker::PhantomData;

pub struct MultilinearKZGProof<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    pub commitment: P::G1,
    pub v: T,
    pub q_taus: Vec<P::G1>,
}

impl<T: PrimeField, P: Pairing> MultilinearKZGProof<T, P> {
    pub fn new(commitment: P::G1, v: T, q_taus: Vec<P::G1>) -> Self {
        Self {
            _marker: PhantomData,
            commitment,
            v,
            q_taus,
        }
    }
}

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
    ) -> MultilinearKZGProof<T, P> {
        let commitment = Self::generate_commitment(polynomial, encrypted_lagrange_basis);
        let opening_points = openings.iter().map(|val| Some(*val)).collect::<Vec<_>>();
        let v_poly = polynomial.evaluate(&opening_points);

        let mut quotients = Vec::with_capacity(openings.len());

        let f_minus_v = polynomial.minus(v_poly.get_evaluation_points().first().unwrap());

        let mut dividend = f_minus_v;

        for (idx, opening) in openings.iter().enumerate() {
            // divide the polynomial by each opening as a factor
            // e.g. if the roots are a = 6, b = 7, c = 0; we divide the polynomial by a - 6, remainder by b - 7 and lastly, c - 0;
            // But in actual fact, we are evaluating the polynomial at the variable points.
            let (mut quotient, remainder) = dividend.compute_quotient_remainder(opening, 0);

            dividend = remainder;
            quotient = MultiLinearPolynomial::blow_up_n_times(
                BlowUpDirection::Left,
                max(idx + 1, openings.len() - (quotient.len().ilog2()) as usize),
                &quotient,
            );

            let quotient_at_tau: P::G1 = Self::evaluate_at_tau(
                &MultiLinearPolynomial::new(&quotient),
                encrypted_lagrange_basis,
            );

            quotients.push(quotient_at_tau);
        }

        MultilinearKZGProof::new(
            commitment,
            *v_poly.get_evaluation_points().first().unwrap(),
            quotients,
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::multilinear::trusted_setup::TrustedSetup;

    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Affine};
    use ark_ec::AffineRepr;

    #[test]
    pub fn test_generate_commitment() {
        let trusted_setup: TrustedSetup<Fr, Bls12_381> =
            TrustedSetup::new(&[Fr::from(5), Fr::from(2), Fr::from(3)]);

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
            &trusted_setup.encrypted_lagrange_basis,
        );

        assert_eq!(
            commitment,
            G1Affine::generator().mul_bigint(Fr::from(42).into_bigint())
        )
    }
}
