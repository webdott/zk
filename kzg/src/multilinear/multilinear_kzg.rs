use ark_ec::pairing::Pairing;
use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;

use ark_ec::PrimeGroup;
use ark_ff::PrimeField;

// Given a set of tau values, we want to generate the lagrange basis array over the boolean hypercube for the number of variables
// We use the check 0 and check one principle
// E.g:
// 000 => (1 - a) * (1 - b) * (1 - c)
// 010 => (1 - a) * (b) * (1 - c)
// 111 => (a) * (b) * (c)
pub fn generate_lagrange_basis_for_n_variables<T: PrimeField>(n: usize, taus: &[T]) -> Vec<T> {
    if n != taus.len() {
        panic!("Length of variables does not match the number of Taus given!")
    }

    let length_of_lagrange_basis = 1 << n;
    let mut lagrange_basis = vec![T::one(); length_of_lagrange_basis];

    for i in 0..length_of_lagrange_basis {
        for (tau_idx, tau) in taus.iter().enumerate() {
            let skip_value = 1 << (n - tau_idx - 1);

            if (i / skip_value) % 2 == 0 {
                // if bit at variable index is 0 (turned off), use check-zero (1 - x)
                lagrange_basis[i] *= T::one() - *tau;
            } else {
                // else use check-one (x)
                lagrange_basis[i] *= tau;
            }
        }
    }

    lagrange_basis
}

// Encrypt lagrange basis by raising the generator to the power of each value in the basis
pub fn encrypt_lagrange_basis<T: PrimeField, P: Pairing>(lagrange_basis: &[T]) -> Vec<P::G1> {
    lagrange_basis
        .iter()
        .map(|&basis_value| P::G1::generator().mul_bigint(basis_value.into_bigint()))
        .collect::<Vec<P::G1>>()
}

// Commitment is gotten by doing an element wise multiplication between encrypted lagrange basis and the multilinear polynomial
pub fn generate_commitment<T: PrimeField, P: Pairing>(
    polynomial: &MultiLinearPolynomial<T>,
    taus: &[T],
) -> P::G1 {
    let no_of_polynomial_variables = polynomial.number_of_variables() as usize;

    if no_of_polynomial_variables != taus.len() {
        panic!("Number of variables does not match the number of Taus given!")
    };

    let lagrange_basis = generate_lagrange_basis_for_n_variables(no_of_polynomial_variables, taus);
    let encrypted_lagrange_basis = encrypt_lagrange_basis::<T, P>(&lagrange_basis);
    let polynomial_evals = polynomial.get_evaluation_points();

    let commitment_array = (0..encrypted_lagrange_basis.len())
        .map(|i| encrypted_lagrange_basis[i].mul_bigint(polynomial_evals[i].into_bigint()))
        .collect::<Vec<_>>();

    commitment_array.iter().sum::<P::G1>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr, G1Affine};
    use ark_ec::AffineRepr;

    #[test]
    pub fn test_lagrange_basis_for_n_variables_with_same_length_of_taus() {
        assert_eq!(
            generate_lagrange_basis_for_n_variables(3, &[Fr::from(5), Fr::from(2), Fr::from(3)]),
            vec![
                Fr::from(-8),
                Fr::from(12),
                Fr::from(16),
                Fr::from(-24),
                Fr::from(10),
                Fr::from(-15),
                Fr::from(-20),
                Fr::from(30),
            ]
        )
    }

    #[test]
    pub fn test_generate_commitment() {
        let commitment = generate_commitment::<Fr, Bls12_381>(
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
            &[Fr::from(5), Fr::from(2), Fr::from(3)],
        );

        assert_eq!(
            commitment,
            G1Affine::generator().mul_bigint(Fr::from(42).into_bigint())
        )
    }
}
